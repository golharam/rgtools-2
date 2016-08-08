package tools;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.InetAddress;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.stream.Collectors;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
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
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/*
 * This program print out coverage information for regions given in a BED file
 * Input: BAM file, BED file of target regions
 * Output: Coverage information for regions in a BED file 
*/
public class CalculateTargetRegionCoverage {
    private static final Log log = Log.getInstance(CalculateTargetRegionCoverage.class);

	public static void main(String[] args) throws IOException {
        if (args.length < 1) {
            System.out.println("Usage: " + CalculateTargetRegionCoverage.class.getCanonicalName() + " bedFile bamFile [outFile]");
            System.exit(1);
        }
        final File bedFile = new File(args[0]);
        final File bamFile = new File(args[1]);
        final File outputFile = args.length >= 3 ? new File(args[2]) : null;
        
        final long start = System.currentTimeMillis();

        log.info("Start with args:" + Arrays.toString(args));
        printConfigurationInfo();

        // open BED file
        BEDCodec bedCodec = new BEDCodec();
        final AbstractFeatureReader<BEDFeature, LineIterator> bedReader = AbstractFeatureReader.getFeatureReader(bedFile.getAbsolutePath(), bedCodec, false);
           
        // open SAM file
        final SamReader samReader = SamReaderFactory.makeDefault().open(bamFile);

       	// open output file
        PrintWriter outWriter = (outputFile != null) ? new PrintWriter(outputFile) : null;
        String outStr = "chr\tstart\tend\tname\tlength\tcoverage\ttotalBases0X\ttotalBases10X";
        if (outWriter != null) {
        	outWriter.println(outStr);
        } else {
        	log.info(outStr);
        }
        
        // iterate BED file
    	long totalReadCount = 0l;
    	long bedRecordCount = 0l;
    	Iterator<BEDFeature> bedIterator = bedReader.iterator();
    	while (bedIterator.hasNext()) {
    		BEDFeature bedFeature = bedIterator.next();

    		// make sure the bedFeature chromosome is in the SAM file header
        	if (samReader.getFileHeader().getSequenceDictionary().getSequence(bedFeature.getContig()) == null) {
        		log.warn("Feature " + bedFeature.getContig() + " does not exist in the SAM reference. Skipping BED feature...");
        		continue;
        	}

    		// iterate of SAM records that overlap this BED feature and record how many SAM reads are each position of the BED feature
        	int bedFeatureStart = bedFeature.getStart();
        	int bedFeatureEnd = bedFeature.getEnd();
        	int bedFeatureLength = bedFeatureEnd - bedFeatureStart + 1;
        	int[] perBaseCoverage = new int[bedFeatureLength];
        	int readCount = 0;
        	
    		SAMRecordIterator samIterator = samReader.query(bedFeature.getContig(), bedFeature.getStart(), bedFeature.getEnd(), false);
    		while (samIterator.hasNext()) {
        		SAMRecord rec = samIterator.next();
        		if (filterRead(rec)) continue;
            	
        		readCount++;
                // walk over each base of the bed feature and add 1 to each base this read covers.
                for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                    final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                    for (int pos=block.getReferenceStart(); pos<=end; ++ pos) {
                        if (pos >= bedFeatureStart && pos <= bedFeatureEnd) {
                        	perBaseCoverage[pos - bedFeatureStart]++;
                        }
                    }
                }
        		
    		}
    		samIterator.close();

    		int totalBases0x = 0;
    		int totalBases10x = 0;
    		double coverage = readCount / (double)bedFeatureLength;
        	for (int i = 0; i < perBaseCoverage.length; i++) {
        		if (perBaseCoverage[i] == 0) {
        			totalBases0x++;
        		}
        		if (perBaseCoverage[i] >= 10) {
        			totalBases10x++;
        		}        		
        	}
        	
        	outStr = String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", bedFeature.getContig(), bedFeature.getStart(), bedFeature.getEnd(), 
        													     bedFeature.getName(), bedFeatureLength, coverage, totalBases0x, totalBases10x);
        	if (outWriter != null)
        		outWriter.println(outStr);
        	else
        		log.info(outStr);

        	totalReadCount += readCount;
    		bedRecordCount++;
        }
    	if (outWriter != null)
    		outWriter.close();
    	samReader.close();
    	bedReader.close();
    	log.info("Found " + totalReadCount + " reads spanning " + bedRecordCount + " BED features");
        final long end = System.currentTimeMillis();
        log.info(String.format("Done. Elapsed time %.3f seconds", (end - start) / 1000.0));                
	}

    /* Return true if the read should be filtered out
     * 
     */
    private static boolean filterRead(SAMRecord rec) {
        // Just plain avoid records that are marked as not-primary
        if (rec.getNotPrimaryAlignmentFlag()) return true;
        // Check for PF reads
        if (rec.getReadFailsVendorQualityCheckFlag()) return true;
        // Check for unmapped reads
        if (rec.getReadUnmappedFlag()) return true;
        // Check for reads that are marked as duplicates
        if (rec.getDuplicateReadFlag()) return true;
        // Don't bother with reads that didn't align uniquely
        if (rec.getMappingQuality() == 0) return true;

        return false;
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
