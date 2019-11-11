package programs;


import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created by chen on 08/02/2015.
 */
public class MahalanobisTest {
    public static void main(String[] args) throws IOException{
        double[][] data = new double[20][3];
        double[][] covarianceMatrix;
        MeshiLineReader reader = new MeshiLineReader("data.dat");
        String line = "";

        int i = 0;
        while ((line = reader.readLine()) != null) {
            String[] words = line.split("  ");
            data[i][0] = Double.valueOf(words[0]) ;
            data[i][1] = Double.valueOf(words[1]) ;
            data[i][2] = Double.valueOf(words[2]) ;
            System.out.println(line+" ; "+words[0]+" , "+words[1]+" , "+words[2]+" "+data[i][0]+" "+data[i][1]+" "+data[i][2]);
            i++;
        }

        double [] mu =  meshi.util.mathTools.Utils.getMeanVector(data) ;
        System.out.println("mu vector");
        for (i = 0; i <mu.length; i++) {
            System.out.print(mu[i]+" ");
        }
        System.out.println();


        covarianceMatrix = meshi.util.mathTools.Utils.getCovarianceMatrix(data);

        System.out.println("Covariance matrix");
        for (i = 0; i < covarianceMatrix.length; i++) {
            for (int j = 0; j < covarianceMatrix.length; j++)
                System.out.print(covarianceMatrix[i][j] + " ");
            System.out.println();
        }

        System.out.println(":) "+meshi.util.mathTools.Utils.Mahalanobis(mu,covarianceMatrix,data[0]));
        double sum = 0;
        for (i = 0; i< data.length; i++) {
            sum += meshi.util.mathTools.Utils.Mahalanobis(mu,covarianceMatrix,data[i]);
        }
        System.out.println(":) :) "+Math.sqrt(sum/data.length));
    }
}
