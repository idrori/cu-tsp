package programs;

import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created by chen on 02/03/2017.
 */
//0123456789012345678901234567890
//ATOM      1  N   GLY A  -2       5.914 -21.810   8.676  1.00 31.05           N
public class MSE2MET {
    public static void main(String[] args) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(args[0]);
        MeshiWriter writer = new MeshiWriter(args[1]);
        String line = "";
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("TER") || line.startsWith("END")) {
                writer.println(line);
                continue;
            }
            if (line.length() < 25) continue;
            if ((line.charAt(24) == '-') ||
                ((line.charAt(24) == ' ') && (line.charAt(25) == '0'))) {
                System.out.println("Removing "+line);
                continue;
            }
            if (line.charAt(13) == 'H') continue;
            if (line.startsWith("ATOM")) writer.println(line);
            else if (line.startsWith("HETATM")) {
                int index = line.indexOf("MSE");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "MET" + line.substring(index + 3);
                    index = line.indexOf(" SE ");
                    if (index < 0) writer.println(newLine);
                    else {
                        newLine = "ATOM  "+line.substring(6,index) + "  SD  MET" + line.substring(index + 9);
                        writer.println(newLine);
                    }
                }
                index = line.indexOf("MLY");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "LYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("MLZ");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "LYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("LLP");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "LYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("OCS");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("CSS");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("LED");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "LEU" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("CME");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("SEC");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("CSO");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("CGV");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("CSD");
                if (index > 0) {
                    String newLine = "ATOM  " + line.substring(6, index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("CSX");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("SNC");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("CAF");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("KST");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "LYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("KCX");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "LYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("K1R");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "CYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("SEP");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "SER" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("TYI");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "TYR" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("DAH");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "PHE" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("TY2");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "TYR" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("LP6");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "LYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
                index = line.indexOf("KPI");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "LYS" + line.substring(index + 3);
                    writer.println(newLine);
                }
            }
        }
        writer.close();
    }
}
