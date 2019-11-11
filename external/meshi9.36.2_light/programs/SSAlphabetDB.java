package programs;

import meshi.energy.EvaluationException;
import meshi.optimizers.OptimizerException;
import meshi.util.MeshiProgram;
import meshi.util.UpdateableException;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.Callable;

/**
 * Created by user on 18/02/2017.
 */
public class SSAlphabetDB implements Callable<Integer> {

    public String[] args= {};
    public static void main(String[] args) throws IOException, OptimizerException, UpdateableException, EvaluationException, Exception {
        String arg2_iniFileName = "-inFileName=";
        String arg3_dssp = "-dsspFile=";
        String arg4_native = "-nativeStructure=";
        String arg5_out = "-outFile=";
        String[] arg = {"-commands=../commands_sidiLocalV2", "-inFileName=QUARK_TS3.pdb", "-dsspFile=QUARK_TS3.dssp", "-nativeStructure=T0755.N.pdb" ,"-outFile=QUARK_TS3.out.pdb", "-seed=1"};

        File folder = new File(".");
        String wDir=folder.getAbsolutePath().substring(0,folder.getAbsolutePath().length()-2);
        System.out.println("Current Working Directory: " + folder.getAbsolutePath());
        File[] listOfFiles = folder.listFiles();
        for (int iDir=0;iDir <listOfFiles.length;iDir++){//targets
            if (listOfFiles[iDir].isDirectory()) {
                System.out.println(listOfFiles[iDir].getName());
                System.setProperty("user.dir", wDir+"\\"+listOfFiles[iDir].getName());

                File target = new File(String.valueOf(wDir+"\\"+listOfFiles[iDir].getName()+"\\pdb"));
                File test = new File(".");
                System.out.println("Current Working Directory: " + test.getAbsolutePath());
                File[] listOfpdbs = target.listFiles();
                for (int iModel=0;iModel <listOfpdbs.length;iModel++) {//targets
                    System.out.println(listOfpdbs[iModel].getName());
                    String commands = "-commands="+wDir+"\\..\\commands_sidiLocalV2";
                    String nwd = test.getAbsolutePath()+"\\";
                    String dsspFilePath=nwd+"dssp\\"+listOfpdbs[iModel].getName().substring(0,listOfpdbs[iModel].getName().lastIndexOf('.'))+".scwrl.dssp";
                    String[] mainArgs = {commands,arg2_iniFileName+nwd+"pdb\\"+listOfpdbs[iModel].getName(),
                            arg3_dssp+dsspFilePath,arg4_native+nwd+"aux_\\"+listOfFiles[iDir].getName()+".N.pdb","-outFile="+nwd+listOfpdbs[iModel].getName()+"out.pdb", "-seed=1"};

                    System.out.println(mainArgs[0]+" "+mainArgs[1]+" "+mainArgs[2]+" "+mainArgs[3]+" "+mainArgs[4]);

                    //MeshiProgram.randomNumberGenerator = null;
                    DsspAlphabetExploration2.main(mainArgs);
                    //Future<Integer> result = null;
                    //ExecutorService executor = Executors.newSingleThreadExecutor();
                    //result = executor.submit(new SSAlphabetDB(mainArgs));
                    //result.get();
                    //DsspAlphabetExploration2.main(mainArgs);
                }
                System.setProperty("user.dir", wDir);
            }
        }
        //DsspAlphabetExploration2.main(arg);
    }

    public SSAlphabetDB(String[] args){
        this.args = args;
    }
    @Override
    public Integer call() {
        try {
            DsspAlphabetExploration2.main(args);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return 1;
    }
}
