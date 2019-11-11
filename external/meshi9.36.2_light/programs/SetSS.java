package programs;

import meshi.util.file.MeshiLineReader;

import java.io.IOException;

public class SetSS {
    public static void main(String[] args)throws IOException{
        MeshiLineReader in = new MeshiLineReader(args[0]);
        int start = Integer.valueOf(args[1]);
        String line;
        int end;
        boolean first = true;
        while((line = in.readLine()) != null) {
                if (!first){
                        System.out.print(",");
                }
                first = false;
                end = start+line.length();
                for (int i = start; i < end ; i++)  {
                        System.out.print((i+1)+"="+line.charAt(i-start));
                        if ((i-start) < line.length()-1)System.out.print(",");
                }
                start = end;
        }
        System.out.println();

    }
}
