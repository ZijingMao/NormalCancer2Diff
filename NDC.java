import java.io.FileInputStream;
import java.io.InputStream;


public class NDC {
	public static String FileNameN;
	public static String FileNameC;
	
	public static void main(String[] args) throws Exception{
		if(args.length != 2){
			System.out.println("Usage: ./Prog FileNameN FileNameC");
			return;
		}
		
		FileNameN = args[0];
		FileNameC = args[1];
		
		System.out.println("Input file: \n"+ FileNameN+"\n"+FileNameC+"\n");
		
		InputStream iN = new FileInputStream(FileNameN);
		InputStream iC = new FileInputStream(FileNameC);
		
		int[][] arrayN = FileHelper.scan(iN, 10);
		int[][] arrayC = FileHelper.scan(iC, 10);
		
		HMM hmm = new HMM();

/*
		double[] PioD = new double[2];
		PioD[0] = 0.9999999999;
		PioD[1] = 1e-10;
		double[][] D2D = new double[2][2];
		D2D[0][0] = 0.9705732139029989;
		D2D[0][1] = 1 - D2D[0][0];
		D2D[1][1] = 0.7143056667193981;
		D2D[1][0] = 1 - D2D[1][1];
		double Weight = 0.35522278;
		int[] path = hmm.decode(arrayC, arrayN, Weight, D2D, PioD);
		
		FileHelper.writeThreeMark(FileHelper.decode(path));

		
		double likelihood = hmm.learn(arrayC, arrayN, 12, 0.001);
		FileHelper.writeMark(hmm.Mark, "mark.txt");
*/

/*
		hmm.DoTwoChains(arrayC, arrayN, 15, 0.001);
		FileHelper.writeMark(hmm.MarkC, "markC.txt");
		FileHelper.writeMark(hmm.MarkN, "markN.txt");
*/

		//FileHelper.writeParam(hmm.PioD, hmm.D2D, hmm.Weight);

		double likelihood = hmm.learn(arrayC, arrayN, 9, 0.001);
		FileHelper.writeGoD(hmm.GoD);

		System.out.println("Successful!");
	}
}
