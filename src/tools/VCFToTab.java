package tools;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VCFToTab extends CommandLineTool {
	private static String version = "16.11.10";
	private static final Log log = Log.getInstance(VCFToTab.class);

	public static void main(String[] args) throws IOException  {
        if (args.length != 2) {
            System.out.println("Usage: " + VCFToTab.class.getCanonicalName() + " vcfFile outFile");
            System.exit(1);
        }
        final File vcfFile = new File(args[0]);
        final File outputFile = new File(args[1]);
        
        final long start = System.currentTimeMillis();

        log.info("Start with args:" + Arrays.toString(args));
        printConfigurationInfo(version);

        // open VCF file
        VCFCodec vcfCodec = new VCFCodec();
        final AbstractFeatureReader<VariantContext, LineIterator> vcfReader = AbstractFeatureReader.getFeatureReader(vcfFile.getAbsolutePath(), vcfCodec, false);

       	// open output file
        PrintWriter outWriter = new PrintWriter(outputFile);

        // Read the header and construct a table of all the possible info fields
        VCFHeader vcfHeader = (VCFHeader) vcfReader.getHeader();

        // Get the format fields
        HashMap<String, VCFFormatHeaderLine> vcfFormatById = new HashMap<String, VCFFormatHeaderLine>();
        HashMap<String, VCFInfoHeaderLine> vcfInfoById = new HashMap<String, VCFInfoHeaderLine>();        
        ArrayList<String> perSampleKeyList = new ArrayList<String>();
        ArrayList<String> keyList = new ArrayList<String>();
        
        // Get format fields
        for (VCFFormatHeaderLine vcfInfo: vcfHeader.getFormatHeaderLines()) {
			vcfFormatById.put(vcfInfo.getID(), vcfInfo);
			perSampleKeyList.add(vcfInfo.getID());
        	
        }

        // Get info fields
        for (VCFInfoHeaderLine vcfInfo : vcfHeader.getInfoHeaderLines()) {
        	vcfInfoById.put(vcfInfo.getID(), vcfInfo);
        	keyList.add(vcfInfo.getID());
        }

        // These are the samples that appear in the vcf file
        ArrayList<String> sampleNames = new ArrayList<String>(vcfHeader.getGenotypeSamples());

        // Print out header fields, followed by info fields, folled by samples.
    	outWriter.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER");
    	for (String key: keyList) {
    		outWriter.write("\t" + key);
    	}
    	// Now print out the genotype format fields for each sample
    	for (String sample : sampleNames) {
    		for (String perSampleKey : perSampleKeyList) {
    			outWriter.write("\t" + sample + "-" + perSampleKey);
    		}
    	}
    	outWriter.write("\n");

    	// print out the variants with the info field values in the same order
        for (final VariantContext vc : vcfReader.iterator()) {
        	String altAlleles = StringUtil.join(",", vc.getAlternateAlleles());
        	String filter = StringUtil.join(",", vc.getFilters());
        	
        	outWriter.write(vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getID() + "\t" + vc.getReference().getDisplayString() + "\t" + 
        			altAlleles + "\t" + vc.getPhredScaledQual() + "\t" + filter);
        	
        	// Print out the info fields
        	for (String key : keyList) {
        		outWriter.write("\t");
        		if (vc.hasAttribute(key)) {
        			List<Object> attrs = vc.getAttributeAsList(key);
        			outWriter.write(StringUtil.join(",", attrs));
        		}
    		}
        	
        	// Print out the genotypes
        	Iterable<Genotype> genotypes = vc.getGenotypesOrderedBy(sampleNames);
        	for (Genotype genotype : genotypes) {
        		for (String perSampleKey : perSampleKeyList) {
        			outWriter.write("\t");
        			if (genotype.hasAnyAttribute(perSampleKey)) {
        				if (perSampleKey.equals("GT")) {
        					outWriter.write(genotype.getGenotypeString());
        				} else if (perSampleKey.equals("PL")) {
        					outWriter.write(genotype.getLikelihoodsString());
        				} else {
        					outWriter.write("" + genotype.getAnyAttribute(perSampleKey));
        				}
        			}
        		}
        	}
        	outWriter.write("\n");
        }

   		outWriter.close();
    	vcfReader.close();
        final long end = System.currentTimeMillis();
        log.info(String.format("Done. Elapsed time %.3f seconds", (end - start) / 1000.0));                
	}

}
