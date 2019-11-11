package meshi.util.mathTools;

import meshi.util.mathTools.Jama.Matrix;

/**
 * Created by chen on 08/02/2015.
 */
public class Utils {
    public static double[][] getCovarianceMatrix(double[][] data) {
        int dim = data[0].length;
        int size = data.length;
        double[][] out = new double[dim][dim];
        double[] mu;
        double[][] diff = new double[size][dim];

        mu = getMeanVector(data);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < size; j++)
                diff[j][i] = data[j][i] - mu[i];
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double sum = 0;
                for (int k = 0; k < size; k++)
                    sum += diff[k][i] * diff[k][j];
                out[i][j] = out[j][i] = sum / (size-1);
            }
        }
        return out;
    }

    public static double[] getMeanVector(double[][] data) {
        int dim = data[0].length;
        int size = data.length;

        double[] mu = new double[dim];
        for (int i = 0; i < dim; i++) {
            double sum = 0;
            for (int j = 0; j < size; j++)
                sum += data[j][i];
            mu[i] = sum / size;
        }
        return mu;
    }

    public static double Mahalanobis(double[] mu, double[][] covariance, double[] datum) {
        int dim = mu.length;
        if (covariance.length != dim) throw new RuntimeException("This is weird");
        if (datum.length != dim) throw new RuntimeException("This is weird");
        for (int i = 0; i < dim; i++)
            if (covariance[i].length != dim) throw new RuntimeException("This is weird");

        Matrix muM = new Matrix(mu);
        Matrix covM = new Matrix(covariance);
        Matrix datumM = new Matrix(datum);
        Matrix diff = datumM.minus(muM);
        Matrix diffT = diff.transpose();
        Matrix out = diff.times(covM.inverse().times(diffT));
        return out.get(0,0);
    }
}
