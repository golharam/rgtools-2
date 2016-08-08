package tools;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.InetAddress;
import java.util.Arrays;
import java.util.Iterator;
import java.util.stream.Collectors;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.FullBEDFeature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

/*
 * This program print out variants from a VCF file for a specified region given in a BED file
 * Input: VCF file, BED file of target regions
 * Output: Variants in target regions 
*/
public class SelectVariants {
    private static final Log log = Log.getInstance(SelectVariants.class);

	public static void main(String[] args) throws IOException {
        if (args.length < 2) {
            System.out.println("Usage: " + SelectVariants.class.getCanonicalName() + " vcfFile bedFile [outFile] [summaryFile] [-exclude dbsnp]");
            System.exit(1);
        }
        final File vcfFile = new File(args[0]);
        final File bedFile = new File(args[1]);
        final File outputFile = args.length >= 3 ? new File(args[2]) : null;
        final File summaryFile = args.length >= 4 ? new File(args[3]) : null;
        
        final long start = System.currentTimeMillis();

        log.info("Start with args:" + Arrays.toString(args));
        printConfigurationInfo();

        // open input VCF file
        VCFCodec vcfCodec = new VCFCodec();
        final AbstractFeatureReader<VariantContext, LineIterator> vcfReader = AbstractFeatureReader.getFeatureReader(vcfFile.getAbsolutePath(), vcfCodec, false);
        loadIndex(vcfFile, vcfCodec);

       	// open output VCF file
        VariantContextWriter vcfWriter = (outputFile != null) ? new VariantContextWriterBuilder().setOutputFile(outputFile).setOutputFileType(VariantContextWriterBuilder.OutputType.VCF).unsetOption(Options.INDEX_ON_THE_FLY).build() : null;
        if (vcfWriter !=  null) {
            vcfWriter.writeHeader((VCFHeader) vcfReader.getHeader());
        }
        
        // open output summary file
        PrintWriter summaryWriter = (summaryFile != null) ? new PrintWriter(summaryFile, "UTF-8") : null;
        if (summaryWriter != null) {
        	summaryWriter.println("#name\tvariant_count");
        }
        
        // open BED file
        AbstractFeatureReader bedReader = AbstractFeatureReader.getFeatureReader(bedFile.getAbsolutePath(), new BEDCodec(), false); // getBEDReader(bedFile);
        
        // TODO: Make sure chromosomes in VCF files and BEd files all match
        
        // now read iterate over the BED file
        final ProgressLogger pl = new ProgressLogger(log, 1000000);
        Iterator<Feature> iter = bedReader.iterator();
        long bedRecordCount = 0l;
        while (iter.hasNext()) {
        	FullBEDFeature bedFeature = (FullBEDFeature) iter.next();
        	log.info("Scanning for variants in region " + bedFeature.getName() + " (" + bedFeature.getContig() + ":" + bedFeature.getStart() + "-" + bedFeature.getEnd() + ")");
        	
            // iterate VCF file
        	long vcfRecordCount = 0l;
        	Iterator<VariantContext> vcIterator = vcfReader.query(bedFeature.getContig(), bedFeature.getStart(), bedFeature.getEnd());
        	while (vcIterator.hasNext()) {
        		VariantContext vc = vcIterator.next();
        		
        		vcfRecordCount++;
            	if (vcfWriter != null) {
                	vcfWriter.add(vc);
                }

            	pl.record(vc.getContig(), vc.getStart());
            }
        	log.info("Found " + vcfRecordCount + " variants");
        	if (summaryWriter != null) {
        		summaryWriter.println(bedFeature.getName() + "\t" + vcfRecordCount);
        	}
            ++bedRecordCount;
        }
        if (vcfWriter != null) {
        	vcfWriter.close();
        }
        if (summaryWriter != null) {
        	summaryWriter.close();
        }																												 																									
        log.info("We saw " + bedRecordCount + " record(s) in file " + bedFile);
        final long end = System.currentTimeMillis();
        log.info(String.format("Done. Elapsed time %.3f seconds", (end - start) / 1000.0));                
	}

    private static void printConfigurationInfo() throws IOException {
        log.info("Executing as " +
                System.getProperty("user.name") + '@' + InetAddress.getLocalHost().getHostName() +
                " on " + System.getProperty("os.name") + ' ' + System.getProperty("os.version") +
                ' ' + System.getProperty("os.arch") + "; " + System.getProperty("java.vm.name") +
                ' ' + System.getProperty("java.runtime.version"));

        log.info(Defaults.allDefaults().entrySet().stream().map(e -> e.getKey() + ':' + e.getValue()).collect(Collectors.<String>joining(" ")));
    }

    private static AbstractFeatureReader getBEDReader(File bedFile) {
	    FeatureCodec bedCodec = new BEDCodec();
	    // get an index
	    Index bedIndex = loadIndex(bedFile, bedCodec);
	    // get a reader
	    AbstractFeatureReader bedReader = AbstractFeatureReader.getFeatureReader(bedFile.getAbsolutePath(), bedCodec, bedIndex);
		return bedReader;
    }
    
    /**
    *
    * @param featureFile the feature file
    * @param codec the codec to decode features with
    * @return an index instance
    */
   public static Index loadIndex(File featureFile, FeatureCodec codec) {
       // lets setup a index file name
       File indexFile = Tribble.indexFile(featureFile);

       // our index instance;
       Index index = null;

       // can we read the index file
       if (indexFile.canRead()) {
    	   log.info("Loading index from disk for index file -> " + indexFile);
           index = IndexFactory.loadIndex(indexFile.getAbsolutePath());
       // else we want to make the index, and write it to disk if possible
       } else {
           log.info("Creating the index and memory, then writing to disk for index file -> " + indexFile);
           index = createAndWriteNewIndex(featureFile,indexFile,codec);
       }

       return index;
   }

   /**
    * creates a new index, given the feature file and the codec
    * @param featureFile the feature file (i.e. .vcf, .bed)
    * @param indexFile the index file; the location we should be writing the index to
    * @param codec the codec to read features with
    * @return an index instance
    */
   public static Index createAndWriteNewIndex(File featureFile, File indexFile, FeatureCodec codec) {
       try {
           Index index = IndexFactory.createLinearIndex(featureFile, codec);

           // try to write it to disk
           LittleEndianOutputStream stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(indexFile)));
           		
           index.write(stream);
           stream.close();

           return index;
       } catch (IOException e) {
           throw new RuntimeIOException("Unable to create index from file " + featureFile,e);
       }
   }

}
