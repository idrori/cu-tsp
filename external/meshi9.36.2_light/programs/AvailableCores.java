package programs;
import java.util.Scanner;

/**
 *  * Created by chen on 09/05/2015.
 *   */
public class AvailableCores {
    static int max32 = 2;
    static int max33 = 2;
    static int max35 = 2;
    static int max223 = 15;

    public static void main(String[] argv) {

        Scanner scanner = new Scanner(System.in);
        int counter = -1;
        String server = null, queue = null;
        String slashString;
        while (scanner.hasNext()) {
            String token = scanner.next();
            if (token.indexOf("sge") > -1) {
                server = token;
                int atIndex = server.indexOf("@");
                queue = server.substring(0, atIndex);
                scanner.next();
                slashString = scanner.next();
            } else continue;
            int slashIndex = slashString.indexOf("/");
            if (slashIndex > 0) {
                String working = slashString.substring(0, slashIndex);
                if (isEnough(server, Integer.parseInt(working))) continue;
                String all = slashString.substring(slashIndex + 1);
                int available = Integer.valueOf(all) - Integer.valueOf(working);
                for (int i = 0; i < available; i++)
                    System.out.println("" + queue);
            }
        }
    }

    public static boolean isEnough(String server, int working) {
        if (server.equals("keasar@sge32")) {
            if (working >= max32) return true;
        }
        if (server.equals("keasar@sge33")) {
            if (working >= max33) return true;
        }
        if (server.equals("keasar@sge35")) {
            if (working >= max35) return true;
        }
        if (server.equals("keasar@sge223")) {
            if (working >= max35) return true;
        }
        Integer i;
        return false;
    }
}