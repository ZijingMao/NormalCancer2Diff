import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;

public class FileHelper {
	
	private static BufferedReader buf = null;
	public final static int MAXLENGTH = 200000; //30804183;
	
	public static int[][] scan(InputStream in, int l) {
		if(in == null){
			System.out.println("No file found.");
			return null;
		}
		
		int[][] array = new int[l][MAXLENGTH];
		
		buf = new BufferedReader(new InputStreamReader(in));
		
		int c, idx = 0, lineNum = 0;
		StringBuilder sb = new StringBuilder();
		do{
			c = read();
			
			if(isDigit(c)){
				sb.append((char)c);
			}else{
				array[lineNum][idx] = Integer.parseInt(sb.toString());
				sb = new StringBuilder();
				if(isLine(c)){
					idx = 0;
					lineNum++;
					while(isLine(peek())){
						read();
					}
				}else if(isComma(c)){
					idx++;
				}
			}
		}while(c != -1 && lineNum < l);
		
		return array;
	}

	private static int peek() {
		int oneByte = -1;
		
		try{
			buf.mark(1);
			oneByte  = buf.read();
			buf.reset();
		}catch(IOException e){
			e.printStackTrace();
		}
		
		return oneByte;
	}

	private static boolean isDigit(int c) {
		if(c >=48 && c < 58){
			return true;
		}
		return false;
	}

	private static int read(){
		int oneByte = -1;
		
		try{
			oneByte  = buf.read();
		}catch(IOException e){
			e.printStackTrace();
		}
		
		return oneByte;
	}
	
	private static boolean isLine(int i){
		if(i =='\n' || i == '\r'){
			return true;
		}
		return false;
	}
	
	private static boolean isComma(int i){
		if(i ==','){
			return true;
		}
		return false;
	}

	public static void writeGoD(double[][] GoD) throws FileNotFoundException {
		OutputStream o = new FileOutputStream("GoD.txt");
		PrintStream ps = new PrintStream(o);
		StringBuilder sb = new StringBuilder();
		for(int i =0; i < GoD.length; i++){
			sb.append(String.format("%.4f", GoD[i][0]));
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()-1);
		ps.println(sb.toString());
		ps.close();
	}

	public static void writeMark(int[] mark, String s) throws FileNotFoundException {
		OutputStream o = new FileOutputStream(s);
		@SuppressWarnings("resource")
		PrintStream ps = new PrintStream(o);
		StringBuilder sb = new StringBuilder();
		for(int i =0; i < mark.length; i++){
			sb.append(mark[i]);
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()-1);
		ps.println(sb.toString());
		ps.close();
	}

	public static void writeParam(double[] pioD, double[][] d2d, double weight) throws FileNotFoundException {
		OutputStream o = new FileOutputStream("params.txt");
		@SuppressWarnings("resource")
		PrintStream ps = new PrintStream(o);
		StringBuilder sb = new StringBuilder();
		
		sb.append("Pi of D: ");
		for(int i =0; i < pioD.length; i++){
			sb.append(pioD[i]);
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()-1);
		sb.append("\n");
		
		sb.append("Transition of D: ");
		for(int i =0; i < d2d.length; i++){
			for(int j = 0; j < d2d[i].length; j++){
				sb.append(pioD[i]);
				sb.append(",");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("\n");
		}
		
		sb.append("Weight: ");
		sb.append(weight);
		sb.append("\n");
		
		ps.println(sb.toString());
		ps.close();
	}

	public static int[][] decode(int[] path) {
		int[][] threeStates = new int[3][path.length];
		
		for(int t = 0; t < path.length; t++){
			threeStates[0][t] = (path[t]&4)/4;
			threeStates[1][t] = (path[t]&2)/2;
			threeStates[2][t] = (path[t]&1)/1;
		}
		
		System.out.println(path.length);
		
		return threeStates;
	}

	public static void writeThreeMark(int[][] decode) throws FileNotFoundException {
		OutputStream o = new FileOutputStream("ThreeMarks.txt");
		@SuppressWarnings("resource")
		PrintStream ps = new PrintStream(o);
		StringBuilder sb;
		for(int i =0; i < decode.length; i++){
			sb = new StringBuilder();
			for(int j = 0; j < MAXLENGTH; j++){
				sb.append(decode[i][j]);
				sb.append(",");
			}
			sb.deleteCharAt(sb.length()-1);
			ps.println(sb.toString());
		}
		ps.close();
	}
	
}
