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

/*
 * This program collected cell and umi metrics for Drop-Seq/Smart-Seq unmapped BAM
 * Input: Cell UMI tagged unmapped BAM file, list of known barcodes
 * Output: # of reads per barcode
*/
public class CollectSingleCellBarcodeMetrics {
    private static String version = "17.08.18";
    private static final Log log = Log.getInstance(CollectSingleCellBarcodeMetrics.class);

    public static void main(String[] args) throws IOException {
        if (args.length < 1) {
            System.out.println("Usage: " + CollectSingleCellBarcodeMetrics.class.getCanonicalName() + " bamFile knownBarcodesFile [outFile]");
            System.exit(1);
        }
        final File bamFile = new File(args[0]);
        final File barcodesFile = new File(args[1]);
        final File outputFile = args.length >= 2 ? new File(args[2]) : null;
        
        final long start = System.currentTimeMillis();

        log.info("Start with args:" + Arrays.toString(args));
        printConfigurationInfo();

        // open SAM file
        final SamReader samReader = SamReaderFactory.makeDefault().open(bamFile);

        // open output file
        PrintWriter outWriter = (outputFile != null) ? new PrintWriter(outputFile) : null;
        
        if (outWriter != null)
            outWriter.close();
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
		log.info("Version " + version);
		String hostname = "";
		try {
			hostname = InetAddress.getLocalHost().getHostName();
		} catch (java.net.UnknownHostException e) {
			hostname = "unknownhost";
		}
	        log.info("Executing as " +
                System.getProperty("user.name") + '@' + hostname +
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
