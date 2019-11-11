package programs;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;

public class getTarget {
    public static final String TARBALL_SUFFIX = ".3D.srv.tar.gz";
    public static void main(String[] argv) throws IOException,URISyntaxException{
        String targetName = argv[0];
        //downloadTargetTarball(targetName);
        downloadTargetNativePdb(targetName);
    }

    public static void downloadTargetTarball(String targetName) throws IOException,URISyntaxException{
        URL url = new URL("http","www.predictioncenter.org","/download_area/CASP8/server_predictions/"+targetName+TARBALL_SUFFIX);
        InputStream inputStream = (InputStream) url.getContent();
        File out = new File(targetName+TARBALL_SUFFIX);
        FileOutputStream outputStream = new FileOutputStream(out);
        int c;
        while ((c = inputStream.read())!= -1) outputStream.write(c);
    }
    public static void downloadTargetNativePdb(String targetName) {

    }


/*
        WebClient webClient = new WebClient();
        HtmlPage page = webClient.getPage(CASP8);
        HtmlElement body = page.getBody();
        List<HtmlTable> tables = body.getHtmlElementsByTagName("table");
        HtmlTable table = tables.get(0);
        List<HtmlTableRow> rows = table.getRows();
        for(HtmlTableRow row : rows) {
            List<HtmlTableCell> cells = row.getCells();
            if(cells.size()> 1)  {
                HtmlTableCell cell = row.getCell(1);
                if (cell.getLocalName().equals("td")) {
                    HtmlAnchor anchor = (HtmlAnchor) cell.getElementsByTagName("a").get(0);
                    if (anchor.getHrefAttribute().equals(targetName+TARBALL_SUFFIX)){
                        System.out.println(anchor);
                        Page p = anchor.openLinkInNewWindow();
                        System.out.println(p);
                    }
                }
            }
        }


//        System.out.println(anchor.getHrefAttribute());

        webClient.closeAllWindows();
  */

}
