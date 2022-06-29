package Benchmarking;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
/*
// readd implementation 'org.jblas:jblas:1.2.5' to gradle
import org.jblas.DoubleMatrix;
*/
/*
import org.ojalgo.matrix.store.PhysicalStore;
import org.ojalgo.matrix.store.Primitive64Store;
*/
/*
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
 */

public class benchmarking {
    // nd4j is slowest however I did not test the implementation with Cuda
    // jblas performance was good because it has a 1D*1D multiply operation, 2D*2D is slow compared to others
    // colt is the fastest due to it having the best suited operations
    // ojalgo has a fast multiply but subtract, transpose, set and fill are incredibly slow compared to everything else
    // I did not test ejml

    public static void main(String[] args){
        double[] mean = new double[]{1.1274211082067633, 2.4200509029594266, 0.10819157839780577, 0.000010494505576886112};
        double[][] covariance = new double[][]{
                {0.5489668520758497, 0.0004185958385229793, -0.009648690178869969, -0.0000008919110937445166},
                {0.00041859583852297785, 0.8982798612543728, 0.05697683281396502, 0.00000593109329120375},
                {-0.009648690178869969, 0.056976832813965025, 0.01345648149181515, 0.0000019015233811744743},
                {-0.0000008919110937445171, 0.00000593109329120375, 0.0000019015233811744743, 0.00000000029405029122465194}
        };
        double[] v = new double[]{1.0525982662842992, 2.497445665993847, 2.031714238116992, 0.3586796851439282};
        double x = 0.0;

        int[] indices = new int[]{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        /*
        DoubleMatrix covariance = new DoubleMatrix(new double[][]{
                {0.5489668520758497, 0.0004185958385229793, -0.009648690178869969, -0.0000008919110937445166},
                {0.00041859583852297785, 0.8982798612543728, 0.05697683281396502, 0.00000593109329120375},
                {-0.009648690178869969, 0.056976832813965025, 0.01345648149181515, 0.0000019015233811744743},
                {-0.0000008919110937445171, 0.00000593109329120375, 0.0000019015233811744743, 0.00000000029405029122465194}
        });
        DoubleMatrix meanMatrix = new DoubleMatrix(new double[] {1.1274211082067633, 2.4200509029594266, 0.10819157839780577, 0.000010494505576886112});
        DoubleMatrix matrix1 = new DoubleMatrix(4, 1);
        DoubleMatrix matrix2 = new DoubleMatrix(4, 4);

        long startTime = System.currentTimeMillis();


        for (int i = 0; i < 100000000; i++) {
            matrix1.put(0, v[0]).put(1, v[1]).put(2, v[2]).put(3, v[3]).subi(meanMatrix);
            double res = matrix2.put(indices, covariance).muliRowVector(matrix1).muliColumnVector(matrix1).sum() * -.5;
            x += res;
        }


        long endTime = System.currentTimeMillis();
        long total = (endTime - startTime) / 1000;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time jblas1: " + total);

         */

        /*
        // add implementation 'org.ojalgo:ojalgo:51.3.0' to gradle build
        // this runs in around 8 seconds, so double colt time
        PhysicalStore.Factory<Double, Primitive64Store> storeFactory = Primitive64Store.FACTORY;
        double[][] data = new double[][]{
                {0.5489668520758497, 0.0004185958385229793, -0.009648690178869969, -0.0000008919110937445166},
                {0.00041859583852297785, 0.8982798612543728, 0.05697683281396502, 0.00000593109329120375},
                {-0.009648690178869969, 0.056976832813965025, 0.01345648149181515, 0.0000019015233811744743},
                {-0.0000008919110937445171, 0.00000593109329120375, 0.0000019015233811744743, 0.00000000029405029122465194}
        };
        Primitive64Store storeCov = storeFactory.rows(data);
        v = new double[]{1.0525982662842992, 2.497445665993847, 2.031714238116992, 0.3586796851439282};
        x = 0.0;
        Primitive64Store storeB = storeFactory.make(1, 4);
        Primitive64Store storeC = storeFactory.make(4, 1);
        Primitive64Store storeD = storeFactory.make(1, 4);
        Primitive64Store storeE = storeFactory.make(1, 1);
        startTime = System.currentTimeMillis();

        for (int i = 0; i < 100000000; i++) {
            // set and fill operations are remarkably slow
            for (int j = 0; j<4; j++){
                storeC.set(j, 0, v[j]-mean[j]);
            }
            storeB.fillRow(0, storeC); // This is faster than transpose since transpose forces a new copy

            storeC.multiply(storeCov, storeD);
            storeD.multiply(storeB, storeE);
            double res = storeE.doubleValue(0) * -0.5;
            x += res;
        }

        endTime = System.currentTimeMillis();
        total = (endTime - startTime) / 1000;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time ojalg1: " + total);
        */

        /*
        matrix1 = new DoubleMatrix(4, 1);
        matrix2 = new DoubleMatrix(4, 4);
        x = 0.0;

        startTime = System.currentTimeMillis();

        for (int i = 0; i < 100000000; i++) {
            matrix1.put(0, v[0]).put(1, v[1]).put(2, v[2]).put(3, v[3]).subi(meanMatrix);
            double res = matrix2.muli(0).addi(covariance).muliRowVector(matrix1).muliColumnVector(matrix1).sum() * -.5;
            x += res;
        }


        endTime = System.currentTimeMillis();
        total = (endTime - startTime) / 1000;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time jblas2: " + total);

         */

        x = 0.0;
        double[][] covarianceA = covariance.clone();
        double[] between = new double[4];

        long startTime = System.currentTimeMillis();

        for (int i = 0; i < 100000000; i++) {
            for (int j = 0; j<v.length; j++){
                between[j] = v[j] - mean[j];
            }

            double t0 = 0;
            for (int j = 0; j<v.length; j++){
//                double t1 = 0;
////                double t1 = between[0] * covarianceA[j][0];
////                t1 += between[1] * covarianceA[j][1];
////                t1 += between[2] * covarianceA[j][2];
////                t1 += between[3] * covarianceA[j][3];
//                for (int k = 0; k<v.length; k++){
//                    t1 += between[k] * covarianceA[j][k];
//                }
//                t0 += between[j] * t1;

                t0 += between[j] * (between[0] * covarianceA[j][0] + between[1] * covarianceA[j][1] + between[2] * covarianceA[j][2] + between[3] * covarianceA[j][3]);
            }
            x += t0 * -0.5;
        }


        long endTime = System.currentTimeMillis();
        long total = (endTime - startTime) / 10;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time with loop0: " + total);

        x = 0.0;
        covarianceA = covariance.clone();
//        double[][] covariancePre = new double[4][4];
//        for (int j = 0; j<v.length; j++){
//            for (int k = 0; k<v.length; k++){
//                covariancePre[j][k] = covarianceA[j][k] * -.5;
//            }
//        }
        between = new double[4];

        startTime = System.currentTimeMillis();

        for (int i = 0; i < 100000000; i++) {
            for (int j = 0; j<v.length; j++){
                between[j] = v[j] - mean[j];
            }

            double t0 = 0;
            for (int j = 0; j<v.length; j++){
                for (int k = 0; k<v.length; k++){
                    t0 += between[k] * between[j] * covarianceA[j][k];
                }
//                t0 += between[j] * t1;
            }
            x += t0 * -0.5;
        }


        endTime = System.currentTimeMillis();
        total = (endTime - startTime) / 10;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time with loop1: " + total);

//        x = 0.0;
//        covarianceA = covariance.toArray2();
//        double[][] covariancePre1 = new double[4][4];
//        double[][] covariancePre2 = new double[4][4];
//        double[][] covariancePre3 = new double[4][4];
//        between = new double[4];
//
//        startTime = System.currentTimeMillis();
//
//        for (int j = 0; j<v.length; j++){
//            for (int k = 0; k<v.length; k++){
////                    t1 += (v[k] - mean[k]) * covarianceA[j][k];
//                covariancePre1[j][k] = mean[k] * covarianceA[j][k];
//                covariancePre2[j][k] = mean[j] * covarianceA[j][k];
//                covariancePre3[j][k] = mean[k] * mean[j] * covarianceA[j][k];
////                t1 += (v[k] * covarianceA[j][k] - mean[k] * covarianceA[j][k]);
//            }
//        }
//
//        for (int i = 0; i < 100000000; i++) {
////            for (int j = 0; j<v.length; j++){
////                between[j] = v[j] - mean[j];
////            }
//
//            double t0 = 0;
//            for (int j = 0; j<v.length; j++){
////                double t1 = 0;
//                for (int k = 0; k<v.length; k++){
////                    t1 += (v[k] - mean[k]) * covarianceA[j][k];
////                    t0 += v[k] * covarianceA[j][k] * v[j] - v[k] * covarianceA[j][k] * mean[j] - covariancePre1[j][k] * v[j] - covariancePre1[j][k] * mean[j];
//                    t0 += v[k] * covarianceA[j][k] * v[j] - v[k] * covariancePre2[j][k] - covariancePre1[j][k] * v[j] - covariancePre3[j][k];
//                }
////                t0 += (v[j] - mean[j]) * t1;
//            }
//            x += t0 * -0.5;
//        }
//
//
//        endTime = System.currentTimeMillis();
//        total = (endTime - startTime) / 10;
//        System.out.println();
//        System.out.println("result: " + x);
//        System.out.println("total time with loop2: " + total);


        x = 0.0;
        covarianceA = covariance.clone();

        double[] covarianceB = new double[16];
        for (int j = 0; j<v.length; j++) {
            System.arraycopy(covarianceA[j], 0, covarianceB, j * v.length, v.length);
        }

        between = new double[4];

        startTime = System.currentTimeMillis();

        for (int i = 0; i < 100000000; i++) {
            for (int j = 0; j<v.length; j++){
                between[j] = v[j] - mean[j];
            }

            double res = 0;
//            for (int j = 0; j<v.length*v.length; j++){
//                res += between[j/v.length] * covarianceB[j] * between[j%v.length];
//            }
            for (int j = 0; j<v.length; j++){
                for (int k = 0; k<v.length; k++){
                    res += covarianceB[j*v.length+k] * between[k] * between[j];
                }
//                res += between[j/v.length] * covarianceB[j] * between[j%v.length];
//                res += covarianceB[j*4] * between[0] * between[j];
//                res += covarianceB[j*4+1] * between[1] * between[j];
//                res += covarianceB[j*4+2] * between[2] * between[j];
//                res += covarianceB[j*4+3] * between[3] * between[j];
            }
            x += res * -0.5;
        }


        endTime = System.currentTimeMillis();
        total = (endTime - startTime) / 10;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time with loop3: " + total);


        DoubleMatrix2D covarianceColt = new DenseDoubleMatrix2D(new double[][]{
                {0.5489668520758497, 0.0004185958385229793, -0.009648690178869969, -0.0000008919110937445166},
                {0.00041859583852297785, 0.8982798612543728, 0.05697683281396502, 0.00000593109329120375},
                {-0.009648690178869969, 0.056976832813965025, 0.01345648149181515, 0.0000019015233811744743},
                {-0.0000008919110937445171, 0.00000593109329120375, 0.0000019015233811744743, 0.00000000029405029122465194}
        });
        DoubleMatrix1D meanMatrixColt1 = new DenseDoubleMatrix1D(new double[]{1.1274211082067633, 2.4200509029594266, 0.10819157839780577, 0.000010494505576886112});

        DoubleMatrix1D matrix1Colt = new DenseDoubleMatrix1D(4);
        DoubleMatrix2D matrix2Colt = new DenseDoubleMatrix2D(4, 4);
        DoubleMatrix1D matrix3Colt = new DenseDoubleMatrix1D(4);
        v = new double[]{1.0525982662842992, 2.497445665993847, 2.031714238116992, 0.3586796851439282};
        x = 0.0;
        Algebra algebra = new Algebra(1.0E-20);

        startTime = System.currentTimeMillis();


        for (int i = 0; i < 100000000; i++) {
            matrix1Colt.assign(v);
            matrix1Colt.assign(meanMatrixColt1, (a, b) -> a-b);

            matrix2Colt.assign(covarianceColt).zMult(matrix1Colt, matrix3Colt);
            double res = algebra.mult(matrix1Colt, matrix3Colt) * -0.5;
            x += res;
        }


        endTime = System.currentTimeMillis();
        total = (endTime - startTime) / 1000;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time colt1: " + total);

        v = new double[]{1.0525982662842992, 2.497445665993847, 2.031714238116992, 0.3586796851439282};
        x = 0.0;

        startTime = System.currentTimeMillis();

        // This is fastest since the assign/set operations are so quick
        for (int i = 0; i < 100000000; i++) {
            matrix1Colt.setQuick(0, v[0]-mean[0]);
            matrix1Colt.setQuick(1, v[1]-mean[1]);
            matrix1Colt.setQuick(2, v[2]-mean[2]);
            matrix1Colt.setQuick(3, v[3]-mean[3]);

            matrix2Colt.assign(covarianceColt).zMult(matrix1Colt, matrix3Colt);
            double res = algebra.mult(matrix1Colt, matrix3Colt) * -0.5;
            x += res;
        }


        endTime = System.currentTimeMillis();
        total = (endTime - startTime) / 1000;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time colt2: " + total);

        /*
        // This is far slower
        matrix1 = new DoubleMatrix(4, 1);
        DoubleMatrix matrix3 = new DoubleMatrix(1, 1);
        DoubleMatrix matrix4 = new DoubleMatrix(1, 4);
        int[] indices2 = new int[]{0,1,2,3};
        x = 0.0;

        startTime = System.currentTimeMillis();

        for (int i = 0; i < 100000000; i++) {
            matrix1.put(0, v[0]).put(1, v[1]).put(2, v[2]).put(3, v[3]).subi(meanMatrix);
            matrix4.put(indices2, matrix1);
            double res = matrix4.mmuli(covariance).mmuli(matrix1, matrix3).get(0) * -.5;
            x += res;
        }


        endTime = System.currentTimeMillis();
        total = (endTime - startTime) / 1000;
        System.out.println();
        System.out.println("result: " + x);
        System.out.println("total time jblas3: " + total);


         */

        /*
        // nd4j is extremely slow without using backend Intel MKL or Cudnn instead of default OpenBLAS
        // and even then not worth it due to JNI overhead for small matrices

        //    Re-add this to gradle build to use nd4j
        //    implementation 'org.nd4j:nd4j-native:1.0.0-M2'
        //    implementation 'org.nd4j:nd4j-api:1.0.0-M2'
        //    implementation 'org.bytedeco.javacpp-presets:openblas:0.3.5-1.4.4'
        //    implementation 'org.bytedeco.javacpp-presets:openblas-platform:0.3.5-1.4.4'
        //    compile "org.nd4j:nd4j-native:1.0.0-M2"
        //    compile group: 'org.deeplearning4j', name: 'deeplearning4j-core', version: '1.0.0-M2'
        //    compile group: 'org.deeplearning4j', name: 'deeplearning4j-nlp', version: '1.0.0-M2'
        //    compile group: 'org.nd4j', name: 'nd4j-native-platform', version: '1.0.0-M2'
        //    compile group: 'org.nd4j', name: 'nd4j-api', version: '1.0.0-M2'
        //    compile 'org.slf4j:slf4j-simple:1.7.36'
        //    compile 'org.slf4j:slf4j-api:1.7.36'

        INDArray covariance1 = Nd4j.create(new double[][]{
                {0.5489668520758497, 0.0004185958385229793, -0.009648690178869969, -0.0000008919110937445166},
                {0.00041859583852297785, 0.8982798612543728, 0.05697683281396502, 0.00000593109329120375},
                {-0.009648690178869969, 0.056976832813965025, 0.01345648149181515, 0.0000019015233811744743},
                {-0.0000008919110937445171, 0.00000593109329120375, 0.0000019015233811744743, 0.00000000029405029122465194}
        });
        INDArray meanMatrix1 = Nd4j.create(new double[] {1.1274211082067633, 2.4200509029594266, 0.10819157839780577, 0.000010494505576886112});
        INDArray matrix11 = Nd4j.create(new double[]{0.0,0.0,0.0,0.0});
        INDArray matrix21 = Nd4j.create(new double[][]{{0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0}});
        v = new double[]{1.0525982662842992, 2.497445665993847, 2.031714238116992, 0.3586796851439282};
        x = 0.0;

        startTime = System.currentTimeMillis();

        for (int i = 0; i < 1; i++) {
            INDArray q = matrix11.putScalar(0, v[0]).putScalar(1, v[1]).putScalar(2, v[2]).putScalar(3, v[3]).subi(meanMatrix1);
            double res = matrix21.muli(0).addi(covariance1).muliRowVector(q).muliColumnVector(q).sum(0,1).getDouble(0);
            x += res;
        }

        endTime = System.currentTimeMillis();
        total = (endTime - startTime) / 1000;
        System.out.println("result: " + x);
        System.out.println();
        System.out.println("total time nd4j: " + total);

         */

        // Did not test ejml
        // according to pre done tests it should be about the same as colt for small matrices
        // http://lessthanoptimal.github.io/Java-Matrix-Benchmark/runtime/2019_02_i53570/
        // implementation 'org.ejml:ejml-all:0.41'
    }
}
