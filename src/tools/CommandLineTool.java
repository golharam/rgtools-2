package tools;

import java.net.InetAddress;
import java.util.stream.Collectors;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.Log;

public class CommandLineTool {
	@SuppressWarnings("unused")
	protected static void printConfigurationInfo(String version) {
		final Log log = Log.getInstance(CommandLineTool.class);

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

}
