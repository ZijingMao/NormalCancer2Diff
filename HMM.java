
public class HMM {

	public double likelihood;
	public int T;
	public int StateoD;
	public double[][][] EoD;
	public double[][] GoD;
	public double[][][][][] fwdoDCNN;
	public double[] scalingC, scalingN, scalingD;   // for scaling forward and backward
	public int StateoCN;
	public double[][] N2X;
	public double[][] N2N;
	public double[] PioN;
	public double[][] C2Y;
	public double[][] C2C;
	public double[] PioC;
	public int SymboloY;
	public int SymboloX;
	public double[][][] PioDCNN;
	public double[][][][][] DCNN;
	public double Weight;
	public double[][] D2D;
	public double[][] CND;
	public double[] PioD;
	public int[] Mark;
	public int[] MarkC;
	public int[] MarkN;
	
	public HMM(){
        StateoD = 2;
        StateoCN = 2;
        SymboloY = 8;
        SymboloX = 8;
        Weight = 0.2;
	 T = FileHelper.MAXLENGTH;

        PioC = new double[StateoCN];
        PioN = new double[StateoCN];
        C2C = new double[StateoCN][StateoCN];
        N2N = new double[StateoCN][StateoCN];
        C2Y = new double[StateoCN][SymboloY];
        N2X = new double[StateoCN][SymboloX];

        PioC[0] = 0.9;
        PioN[0] = 0.9;
        for (int i = 1; i < StateoCN; i++)
        {
            PioC[i] = 0.1 / (StateoCN - 1);
            PioN[i] = 0.1 / (StateoCN - 1);
        }

        for (int i = 0; i < StateoCN; i++)
        {
            for (int j = 0; j < StateoCN; j++)
            {
                if (i == j)
                {
                    C2C[i][j] = 0.90;
                    N2N[i][j] = 0.90;
                }
                else
                {
                    C2C[i][j] = 0.10 / (StateoCN - 1);
                    N2N[i][j] = 0.10 / (StateoCN - 1);
                }
            }
        }

        for (int i = 0; i < StateoCN; i++)
        {
            for (int j = 0; j < SymboloX; j++)
            {
                C2Y[i][j] = 1.0 / SymboloX;
                N2X[i][j] = 1.0 / SymboloX;
            }
        }

        PioD = new double[StateoD];
        D2D = new double[StateoD][StateoD];

        PioD[0] = 0.9;
        for (int i = 1; i < StateoD; i++)
        {
            PioD[i] = 0.1 / (StateoD - 1);
        }

        for (int i = 0; i < StateoD; i++)
        {
            for (int j = 0; j < StateoD; j++)
            {
                if (i == j)
                {
                    D2D[i][j] = 0.7;
                }
                else
                {
                    D2D[i][j] = 0.3 / (StateoD - 1);
                }
            }
        }

        CND = new double[StateoD][StateoCN * 2];
        CND[0][0] = 1;
        CND[0][3] = 1;
        CND[0][1] = 0;
        CND[0][2] = 0;
        CND[1][0] = 0;
        CND[1][3] = 0;
        CND[1][1] = 1;
        CND[1][2] = 1;
    }
	
	/*
	public int[] decode(int[][] arrayC, int[][] arrayN, double Weight, double[][] D2D, double[] PioD) throws Exception{
		if (arrayC == null || arrayN == null)
            throw new Exception("Could not load observed data.");

        if (arrayC.length == 0 || arrayN.length == 0) {
            likelihood = 0.0;
            return null;
        }
        
        T = FileHelper.MAXLENGTH;
        
        DoTwoChains(arrayC, arrayN, 12, 0.001);
        
        MarkC = decodeC(arrayC);
        
        MarkN = decodeN(arrayN);
        
        // Viterbi-forward algorithm.
        int minState;
        double minWeight;
        double weight;

        int[][] s = new int[T][StateoD];
        double[][] a = new double[T][StateoD];
        
        for (int i = 0; i < StateoD; i++) {
    		double p = 0.0;
    		for (int j_c = 0; j_c < arrayC.length; j_c++){
                p += Math.log(C2Y[MarkC[0]][arrayC[j_c][0]]);
            }
    		for (int j_n = 0; j_n < arrayN.length; j_n++){
                p += Math.log(N2X[MarkN[0]][arrayN[j_n][0]]);
            }
            		
    		a[0][i] = (- Math.log(PioD[i])
                    - Math.log(Weight * PioN[MarkN[0]] + (1 - Weight) * CND[MarkN[0]][MarkC[0]+i*2])
                    - Math.log(PioC[MarkC[0]]) - p);
        }
        
        // Induction
        for(int t = 1; t < T; t++){
        	for(int j = 0; j < StateoD; j++){
        		minState = 0;
        		minWeight = a[t-1][0] - Math.log(D2D[0][j]) - Math.log(C2C[0][MarkC[t]])
        				- Math.log(Weight * N2N[0][MarkN[t]] + (1 - Weight) * CND[MarkN[t]][MarkC[t]+j*2]);
		
        		for(int i = 0; i < StateoD; i++){
        			weight = a[t-1][i] - Math.log(D2D[i][j]) - Math.log(C2C[MarkC[t-1]][MarkC[t]])
        					- Math.log(Weight * N2N[MarkN[t-1]][MarkN[t]] + (1 - Weight) * CND[MarkN[t]][MarkC[t]+j*2]);
        			if(weight < minWeight){
        				minState = i;
            			minWeight = weight;
        			}
        		}
		
        		s[t][j] = minState;
        		double p = 0.0;
        		for (int j_c = 0; j_c < arrayC.length; j_c++){
                    p += Math.log(C2Y[MarkC[t]][arrayC[j_c][t]]);
                }
        		for (int j_n = 0; j_n < arrayN.length; j_n++){
                     += Math.log(N2X[MarkN[t]][arrayN[j_n][t]]);
                }
        		
        		a[t][j] = (minWeight - p);
        	}
        }
        
        // Find minimum value for time T-1
        minState = 0;
        minWeight = a[T-1][0];
        
        for(int j = 0; j < StateoD; j++){
    		if (a[T - 1][j] < minWeight)
            {
                minState = j;
                minWeight = a[T - 1][j];
            }
        }
        
        // Trackback
        int[] path = new int[T];
        path[T - 1] = minState;
        
        for (int t = T - 2; t >= 0; t--)
            path[t] = s[t + 1][path[t + 1]];
        
        // Returns the sequence probability as an out parameter
        likelihood = -minWeight;
        
        // Returns the most likely (Viterbi path) for the given sequence
        return path;
	}
	*/
	
	
	public int[] decode(int[][] arrayC, int[][] arrayN, double Weight, double[][] D2D, double[] PioD) throws Exception{
		if (arrayC == null || arrayN == null)
            throw new Exception("Could not load observed data.");

        if (arrayC.length == 0 || arrayN.length == 0) {
            likelihood = 0.0;
            return new int[0];
        }
        
        T = FileHelper.MAXLENGTH;
        
        DoTwoChains(arrayC, arrayN, 20, 0.001);
        
        // Viterbi-forward algorithm.
        int minStateoD, minStateoC, minStateoN;
        double minWeight;
        double weight;

        int[][][][] s = new int[T][StateoD][StateoCN][StateoCN];
        double[][][][] a = new double[T][StateoD][StateoCN][StateoCN];
        
        // Base
        for (int i = 0; i < StateoD; i++) {
            for(int i_c = 0; i_c < StateoCN; i_c++){
            	for(int i_n = 0; i_n < StateoCN; i_n++){
            		double p = 0.0;
            		for (int j_c = 0; j_c < arrayC.length; j_c++){
                        p += Math.log(C2Y[i_c][arrayC[j_c][0]]);
                    }
            		for (int j_n = 0; j_n < arrayN.length; j_n++){
                        p += Math.log(N2X[i_n][arrayN[j_n][0]]);
                    }
            		
            		a[0][i][i_c][i_n] = (- Math.log(PioD[i])
                            - Math.log(Weight * PioN[i_n] + (1 - Weight) * CND[i_n][i_c+i*2])
                            - Math.log(PioC[i_c]) - p);
            	}
            }
        }
        
        // Induction
        for(int t = 1; t < T; t++){
        	for(int j = 0; j < StateoD; j++){
        		for(int l_c = 0; l_c < StateoCN; l_c++){
                	for(int l_n = 0; l_n < StateoCN; l_n++){
                		minStateoD = 0;
                		minStateoC = 0;
                		minStateoN = 0;
                		minWeight = a[t-1][0][0][0] - Math.log(D2D[0][j]) - Math.log(C2C[0][l_c])
                				- Math.log(Weight * N2N[0][l_n] + (1 - Weight) * CND[l_n][l_c+j*2]);
        		
		        		for(int i = 0; i < StateoD; i++){
		        			for(int i_c = 0; i_c < StateoCN; i_c++){
		                    	for(int i_n = 0; i_n < StateoCN; i_n++){
				        			weight = a[t-1][i][i_c][i_n] - Math.log(D2D[i][j]) - Math.log(C2C[i_c][l_c])
				        					- Math.log(Weight * N2N[i_n][l_n] + (1 - Weight) * CND[l_n][l_c+j*2]);
				        			if(weight < minWeight){
				        				minStateoD = i;
				        				minStateoC = i_c;
				        				minStateoN = i_n;
				            			minWeight = weight;
				        			}
		                    	}
		                    }
		        		}
        		
		        		s[t][j][l_c][l_n] = minStateoD*4 + minStateoC*2 + minStateoN;
		        		double p = 0.0;
                		for (int j_c = 0; j_c < arrayC.length; j_c++){
                            p += Math.log(C2Y[l_c][arrayC[j_c][t]]);
                        }
                		for (int j_n = 0; j_n < arrayN.length; j_n++){
                            p += Math.log(N2X[l_n][arrayN[j_n][t]]);
                        }
                		
                		a[t][j][l_c][l_n] = (minWeight - p);
                	}
        		}
        	}
        }
        
        // Find minimum value for time T-1
        minStateoD = 0;
        minStateoC = 0;
        minStateoN = 0;
        minWeight = a[T-1][0][0][0];
        
        for(int j = 0; j < StateoD; j++){
    		for(int l_c = 0; l_c < StateoCN; l_c++){
            	for(int l_n = 0; l_n < StateoCN; l_n++){
            		if (a[T - 1][j][l_c][l_n] < minWeight)
                    {
                        minStateoD = j;
                        minStateoC = l_c;
                        minStateoN = l_n;
                        minWeight = a[T - 1][j][l_c][l_n];
                    }
            	}
    		}
        }
        
        // Trackback
        int[] path = new int[T];
        path[T - 1] = minStateoD*4 + minStateoC*2 + minStateoN;
        
        for (int t = T - 2; t >= 0; t--)
            path[t] = s[t + 1][(path[t + 1]&4)/4][(path[t + 1]&2)/2][(path[t + 1]&1)/1];
        
        // Returns the sequence probability as an out parameter
        likelihood = -minWeight;
        
        // Returns the most likely (Viterbi path) for the given sequence
        return path;
	}

	
	private int[] decodeN(int[][] observations) {
		// Viterbi-forward algorithm.
		T = FileHelper.MAXLENGTH;
        int states = StateoCN;
        int minState;
        double minWeight;
        double weight;

        double[] pi = PioN;
        double[][] A = N2N;

        int[][] s = new int[StateoCN][T];
        double[][] a = new double[StateoCN][T];

        // Base
        for (int i = 0; i < states; i++)
        {
            a[i][0] = -Math.log(pi[i]);
            for (int j = 0; j < observations.length; j++)
                a[i][0] -= Math.log(N2X[i][observations[j][0]]);
        }

        // Induction
        for (int t = 1; t < T; t++)
        {
            int observation = 0;
            for (int j = 0; j < states; j++)
            {
                minState = 0;
                minWeight = a[0][t - 1] - Math.log(A[0][j]);

                for (int i = 1; i < states; i++)
                {
                    weight = a[i][t - 1] - Math.log(A[i][j]);

                    if (weight < minWeight)
                    {
                        minState = i;
                        minWeight = weight;
                    }

                }
                a[j][t] = minWeight;
                for (int k = 0; k < observations.length; k++)
                {
                    observation = observations[k][t];
                    a[j][t] = a[j][t] - Math.log(N2X[j][observation]);
                }
                s[j][t] = minState;
            }
        }

        // Find minimum value for time T-1
        minState = 0;
        minWeight = a[0][T - 1];

        for (int i = 1; i < states; i++)
        {
            if (a[i][T - 1] < minWeight)
            {
                minState = i;
                minWeight = a[i][T - 1];
            }
        }


        // Trackback
        int[] path = new int[T];
        
        path[T - 1] = minState;

        for (int t = T - 2; t >= 0; t--)
        {
            path[t] = s[path[t + 1]][t + 1];
        }
        
        // Returns the most likely (Viterbi path) for the given sequence
        return path;
	}

	private int[] decodeC(int[][] observations) {
		// Viterbi-forward algorithm.
		T = FileHelper.MAXLENGTH;
        int states = StateoCN;
        int minState;
        double minWeight;
        double weight;

        double[] pi = PioC;
        double[][] A = C2C;

        int[][] s = new int[StateoCN][T];
        double[][] a = new double[StateoCN][T];

        // Base
        for (int i = 0; i < states; i++)
        {
            a[i][0] = -Math.log(pi[i]);
            for (int j = 0; j < observations.length; j++)
                a[i][0] -= Math.log(C2Y[i][observations[j][0]]);
        }

        // Induction
        for (int t = 1; t < T; t++)
        {
            int observation = 0;
            for (int j = 0; j < states; j++)
            {
                minState = 0;
                minWeight = a[0][t - 1] - Math.log(A[0][j]);

                for (int i = 1; i < states; i++)
                {
                    weight = a[i][t - 1] - Math.log(A[i][j]);

                    if (weight < minWeight)
                    {
                        minState = i;
                        minWeight = weight;
                    }

                }
                a[j][t] = minWeight;
                for (int k = 0; k < observations.length; k++)
                {
                    observation = observations[k][t];
                    a[j][t] = a[j][t] - Math.log(C2Y[j][observation]);
                }
                s[j][t] = minState;
            }
        }

        // Find minimum value for time T-1
        minState = 0;
        minWeight = a[0][T - 1];

        for (int i = 1; i < states; i++)
        {
            if (a[i][T - 1] < minWeight)
            {
                minState = i;
                minWeight = a[i][T - 1];
            }
        }


        // Trackback
        int[] path = new int[T];
        
        path[T - 1] = minState;

        for (int t = T - 2; t >= 0; t--)
        {
            path[t] = s[path[t + 1]][t + 1];
        }
        
        // Returns the most likely (Viterbi path) for the given sequence
        return path;
	}

	public double learn(int[][] arrayC, int[][] arrayN, int iteration, double tolerance)
		throws Exception{
		if(tolerance == 0 || iteration == 0){
			throw new Exception("Iterations, tolerance should not be 0.");
		}
		
		EoD = new double[T][StateoD][StateoD];
        GoD = new double[T][StateoD];
        
        int currentIteration = 1;
        boolean stop = false;
        
        // old likelihood and new likelihood in order to check converge
        double oldLikelihoodD = Double.MIN_VALUE;
        double newLikelihoodD = 0;
        
        DoTwoChains(arrayC, arrayN, 12, tolerance);

	FileHelper.writeMark(MarkC, "markC.txt");
	FileHelper.writeMark(MarkN, "markN.txt");
        
        PioDCNN = new double[StateoD][StateoCN][StateoCN];
        DCNN = new double[T][StateoD][StateoCN][StateoCN][StateoCN];
        
        // do iteration for CND layer
        System.out.println("Start estimate differential...");
        
        do // Until convergence or max iterations
        {
            // For each sequence in the array input
            System.out.println("Current iteration: "+currentIteration);

            double[][][][] fwdoD = forwardD(arrayC, arrayN);
            double[][][][] bwdoD = backwardD(arrayC, arrayN);

            System.out.println("Calculate gamma and epsilon values for next computations.");
            for (int t = 0; t < T; t++)
            {
                double s = 0.0;                    
                for (int k = 0; k < StateoD; k++)
                {
                    double sum = 0.0;
                    for (int i = 0; i < StateoCN; i++)
                    {
                        for (int j = 0; j < StateoCN; j++)
                        {
                            sum += Math.exp(fwdoD[t][k][i][j] + bwdoD[t][k][i][j]);
                        }
                    }
                    s += GoD[t][k] = sum;
                }

                if (s != 0) // Scaling for computing
                {
                    for (int k = 0; k < StateoD; k++)
                    {
                        GoD[t][k] /= s;
                    }
                }
            }

            
            double m = 0;
            for (int k = 0; k < StateoD; k++)
            {
                for (int i_c = 0; i_c < StateoCN; i_c++)
                {
                    for (int i_n = 0; i_n < StateoCN; i_n++)
                    {
                        m += PioDCNN[k][i_c][i_n] = Math.exp(fwdoD[0][k][i_c][i_n] + bwdoD[0][k][i_c][i_n]);
                    }
                }
            }

            for (int k = 0; k < StateoD; k++)
            {
                for (int i_c = 0; i_c < StateoCN; i_c++)
                {
                    for (int i_n = 0; i_n < StateoCN; i_n++)
                    {
                        PioDCNN[k][i_c][i_n] /= m;
                    }
                }
            }

            for (int t = 1; t < T; t++)
            {
                double s = 0.0;
                for (int k = 0; k < StateoD; k++)
                {
                    for (int i_c = 0; i_c < StateoCN; i_c++)
                    {
                        for (int i_n = 0; i_n < StateoCN; i_n++)
                        {
                            for (int l_n = 0; l_n < StateoCN; l_n++)
                            {
                                s += DCNN[t][k][i_c][i_n][l_n] = Math.exp(fwdoDCNN[t][k][i_c][i_n][l_n] + bwdoD[t][k][i_c][i_n]);
                            }
                        }
                    }
                }

                if (s != 0) // Scaling for computing
                {
                    for (int k = 0; k < StateoD; k++)
                    {
                        for (int i_c = 0; i_c < StateoCN; i_c++)
                        {
                            for (int i_n = 0; i_n < StateoCN; i_n++)
                            {
                                for (int l_n = 0; l_n < StateoCN; l_n++)
                                {
                                    DCNN[t][k][i_c][i_n][l_n] /= s;
                                }
                            }
                        }
                    }
                }
            }

            // Calculate epsilon values for next computations
            for (int t = 0; t < T - 1; t++)
            {
                double s = 0.0;
                for (int k = 0; k < StateoD; k++)
                {
                    for (int l = 0; l < StateoD; l++)
                    {
                        double sum = 0.0;
                        for (int i_c = 0; i_c < StateoCN; i_c++)
                        {
                            for (int i_n = 0; i_n < StateoCN; i_n++)
                            {
                                double p = 0.0;
                                for (int j_c = 0; j_c < arrayC.length; j_c++)
                                {
                                    p += Math.log(C2Y[i_c][arrayC[j_c][t + 1]]);
                                }
                                for (int j_n = 0; j_n < arrayN.length; j_n++)
                                {
                                    p += Math.log(N2X[i_n][arrayN[j_n][t + 1]]);
                                }
                                for (int l_c = 0; l_c < StateoCN; l_c++)
                                {
                                    for (int l_n = 0; l_n < StateoCN; l_n++)
                                    {
                                        sum += Math.exp(fwdoD[t][k][l_c][l_n]
                                            + Math.log(D2D[k][l]) + bwdoD[t + 1][l][i_c][i_n] + Math.log(C2C[l_c][i_c])
                                            + Math.log(Weight * N2N[l_n][i_n] + (1 - Weight) * CND[i_n][i_c + l * 2]) + p);
                                    }
                                }
                            }
                        }
                        EoD[t][k][l] = sum;
                        s += EoD[t][k][l];
                    }
                }

                // scaling
                if (s != 0)
                {
                    for (int k = 0; k < StateoD; k++)
                    {
                        for (int l = 0; l < StateoD; l++)
                        {
                            EoD[t][k][l] /= s;
                        }
                    }
                }                   
            }

            // Compute log-likelihood for the given sequence
            for (int t = 0; t < T; t++)
                newLikelihoodD += scalingD[t];

            System.out.println("Value:");
            System.out.println(Weight);
            System.out.println(newLikelihoodD);
            System.out.println(D2D[0][0]);
            System.out.println(D2D[1][1]);

            if (checkConvergence(oldLikelihoodD, newLikelihoodD,
                currentIteration, iteration, tolerance))
            {
                stop = true;
            }
            else
            {
                // Continue with parameter re-estimation
                currentIteration++;
                oldLikelihoodD = newLikelihoodD;
                newLikelihoodD = 0.0;

                // Re-estimation of initial state probabilities
                for (int k = 0; k < StateoD; k++)
                {
                    PioD[k] = GoD[0][k];
                    if (PioD[k] < 1e-10)
                    {
                        PioD[k] = 1e-10;
                    }
                    if (PioD[k] > 1 - 1e-10)
                    {
                        PioD[k] = 1 - 1e-10;
                    }
                }

                // Re-estimation of transition probabilities 
                for (int i = 0; i < StateoD; i++)
                {
                    for (int j = 0; j < StateoD; j++)
                    {
                        double den = 0, num = 0;

                        for (int l = 0; l < T - 1; l++)
                            num += EoD[l][i][j];
                        for (int l = 0; l < T - 1; l++)
                            den += GoD[l][i];

                        D2D[i][j] = (den != 0) ? num / den : 0.0;
                    }
                }

                System.out.println("Re-estimation of Weight...");
                //Weight = 0.3;
                Weight = WeightEstimation(arrayC, arrayN);
            }
        } while (!stop);

        for (int i = 0; i < 1; i++)
        {
            this.Mark = markState(GoD);
        }

        return newLikelihoodD;
	}

	public int[] markState(double[][] gamma) {
		 // Get the length of the sequence
        int len = T;

        int[] mark = new int[len];

        // Store the mark of sequence using max likelihood
        for (int i = 0; i < len; i++)
        {
            mark[i] = 0;
            for (int j = 0; j < StateoD - 1; j++)
            {
                if (gamma[i][mark[i]] < gamma[i][j + 1])
                {
                    mark[i] = j + 1;
                }
            }
        }

        return mark;
	}

	public double WeightEstimation(int[][] arrayC, int[][] arrayN) {
		double[][][][] Denominator = new double[StateoD][StateoCN][StateoCN][StateoCN];
        double[][][] D0 = new double[StateoD][StateoCN][StateoCN];
        double[][][][] Numerator = new double[StateoD][StateoCN][StateoCN][StateoCN];
        double[][][] N0 = new double[StateoD][StateoCN][StateoCN];
        for (int k = 0; k < StateoD; k++)
        {
            for (int i_c = 0; i_c < StateoCN; i_c++)
            {
                for (int i_n = 0; i_n < StateoCN; i_n++)
                {
                    for (int l_n = 0; l_n < StateoCN; l_n++)
                    {
                        Denominator[k][i_c][i_n][l_n] = CND[i_n][i_c + k * 2] / (N2N[l_n][i_n] - CND[i_n][i_c + k * 2]);
                    }
                    D0[k][i_c][i_n] = CND[i_n][i_c + k * 2] / (PioN[i_n] - CND[i_n][i_c + k * 2]);
                }
            }
        }

        for (int k = 0; k < StateoD; k++)
        {
            for (int i_c = 0; i_c < StateoCN; i_c++)
            {
                for (int i_n = 0; i_n < StateoCN; i_n++)
                {
                    N0[k][i_c][i_n] = PioDCNN[k][i_c][i_n];
                }
            }
        }
       
        for (int t = 1; t < T; t++)
        {
            for (int k = 0; k < StateoD; k++)
            {
                for (int i_c = 0; i_c < StateoCN; i_c++)
                {
                    for (int i_n = 0; i_n < StateoCN; i_n++)
                    {
                        for (int l_n = 0; l_n < StateoCN; l_n++)
                        {
                            Numerator[k][i_c][i_n][l_n] += DCNN[t][k][i_c][i_n][l_n];
                        }
                    }
                }
            }
        }

        double f = F(0.0001, N0, D0, Numerator, Denominator);
        double w = 0.009;

        for (int i = 0; i < 100; i++)
        {
            if (f * F(w, N0, D0, Numerator, Denominator) <= 0)
            {
                f = F(w, N0, D0, Numerator, Denominator);
                for (int j = 0; j < 100; j++)
                {
                    w -= 0.0001;
                    if (f * F(w, N0, D0, Numerator, Denominator) <= 0)
                    {
                        f = F(w, N0, D0, Numerator, Denominator);
                        for (int k = 0; k < 100; k++)
                        {
                            w += 0.000001;
                            if (f * F(w, N0, D0, Numerator, Denominator) <= 0)
                            {
                                f = F(w, N0, D0, Numerator, Denominator);
                                for (int m = 0; m < 100; m++)
                                {
                                    w -= 0.00000001;
                                    if (f * F(w, N0, D0, Numerator, Denominator) <= 0)
                                    {
                                        return w;
                                    }
                                    f = F(w, N0, D0, Numerator, Denominator);
                                }
                            }
                            f = F(w, N0, D0, Numerator, Denominator);
                        }
                    }
                    f = F(w, N0, D0, Numerator, Denominator);
                }
            }
            f = F(w, N0, D0, Numerator, Denominator);
            if (i < 99)
            {
                w += 0.01;
            }
            else
            {
                w = 1e-10;
            }
        }

        return w;
	}

	public double F(double w, double[][][] N0, double[][][] D0, double[][][][] Numerator, double[][][][] Denominator)
    {
        double v = 0;

        for (int k = 0; k < StateoD; k++)
        {
            for (int i_c = 0; i_c < StateoCN; i_c++)
            {
                for (int i_n = 0; i_n < StateoCN; i_n++)
                {
                    v += N0[k][i_c][i_n] / (w + D0[k][i_c][i_n]);
                }
            }
        }

        for (int k = 0; k < StateoD; k++)
        {
            for (int i_c = 0; i_c < StateoCN; i_c++)
            {
                for (int i_n = 0; i_n < StateoCN; i_n++)
                {
                    for (int l_n = 0; l_n < StateoCN; l_n++)
                    {
                        v += Numerator[k][i_c][i_n][l_n] / (w + Denominator[k][i_c][i_n][l_n]);
                    }
                }
            }
        }

        return v;
    }
	
	public double[][][][] backwardD(int[][] arrayC, int[][] arrayN) {
		double[][][][] bwd = new double[T][StateoD][StateoCN][StateoCN];

        for (int k = 0; k < StateoD; k++)
        {
            for (int i_c = 0; i_c < StateoCN; i_c++)
            {
                for (int i_n = 0; i_n < StateoCN; i_n++)
                {
                    bwd[T - 1][k][i_c][i_n] = -scalingD[T - 1];
                }
            }
        }

        System.out.println("Backward...");

        for (int t = T - 2; t >= 0; t--)
        {
            for (int k = 0; k < StateoD; k++)
            {
                for (int i_c = 0; i_c < StateoCN; i_c++)
                {
                    for (int i_n = 0; i_n < StateoCN; i_n++)
                    {
                        double sum = 0.0;
                        for (int l_d = 0; l_d < StateoD; l_d++)
                        {
                            for (int l_c = 0; l_c < StateoCN; l_c++)
                            {
                                for (int l_n = 0; l_n < StateoCN; l_n++)
                                {
                                    double p = 0.0;
                                    for (int j_c = 0; j_c < arrayC.length; j_c++)
                                    {
                                        p += Math.log(C2Y[l_c][arrayC[j_c][t + 1]]);
                                    }
                                    for (int j_n = 0; j_n < arrayN.length; j_n++)
                                    {
                                        p += Math.log(N2X[l_n][arrayN[j_n][t + 1]]);
                                    }
                                    sum += Math.exp((bwd[t + 1][l_d][l_c][l_n] + Math.log(D2D[k][l_d]) + Math.log(C2C[i_c][l_c])
                                        + Math.log(Weight * N2N[i_n][l_n] + (1 - Weight) * CND[l_n][l_c + l_d * 2]) + p));
                                }
                            }
                        }
                        bwd[t][k][i_c][i_n] += Math.log(sum) - scalingD[t];
                    }
                }
            }
        }

        return bwd;
	}

	public double[][][][] forwardD(int[][] arrayC, int[][] arrayN) {
		double[][][][] fwd = new double[T][StateoD][StateoCN][StateoCN];

        fwdoDCNN = new double[T][StateoD][StateoCN][StateoCN][StateoCN];
        scalingD = new double[T];

        // 1. Initialization
        for (int k = 0; k < StateoD; k++)
        {
            for (int i_c = 0; i_c < StateoCN; i_c++)
            {
                for (int i_n = 0; i_n < StateoCN; i_n++)
                {
                    double p = 0.0;
                    
                    
                    for (int j_n = 0; j_n < arrayN.length; j_n++)
                    {
                        p += Math.log(N2X[i_n][arrayN[j_n][0]]);
                    }
                    for (int j_c = 0; j_c < arrayC.length; j_c++)
                    {
                        p += Math.log(C2Y[i_c][arrayC[j_c][0]]);
                    }
                    fwd[0][k][i_c][i_n] = Math.log(PioD[k]) + Math.log(PioC[i_c]) +
                        Math.log(Weight * PioN[i_n] + (1 - Weight) * CND[i_n][i_c + k * 2]) + p;

                    //scalingD[0] += fwd[0][k][i_c][i_n];
                }
            }
        }

        double y = 0;

        for (int k = 0; k < StateoD; k++)
        {
            for (int i_c = 0; i_c < StateoCN; i_c++)
            {
                for (int i_n = 0; i_n < StateoCN; i_n++)
                {
                	y += Math.exp(-fwd[0][0][0][0]+fwd[0][k][i_c][i_n]);
                }
            }
        }


        scalingD[0] = Math.log(y) + fwd[0][0][0][0];


        // scaling
        if (scalingD[0] != 0)
        {
            for (int k = 0; k < StateoD; k++)
            {
                for (int i_c = 0; i_c < StateoCN; i_c++)
                {
                	
                    for (int i_n = 0; i_n < StateoCN; i_n++)
                    {
                        fwd[0][k][i_c][i_n] = fwd[0][k][i_c][i_n] - scalingD[0];
                    }
                }
            }
        }

        System.out.println("Forward...");
         // 2. Induction
        for (int t = 1; t < T; t++)
        {
            for (int k = 0; k < StateoD; k++)
            {
                for (int i_c = 0; i_c < StateoCN; i_c++)
                {
                    for (int i_n = 0; i_n < StateoCN; i_n++)
                    {
                        double p = 0.0;
                        for (int j_n = 0; j_n < arrayN.length; j_n++)
                        {
                            p += Math.log(N2X[i_n][arrayN[j_n][t]]);
                        }
                        for (int j_c = 0; j_c < arrayC.length; j_c++)
                        {
                            p += Math.log(C2Y[i_c][arrayC[j_c][t]]);
                        }

                        double[] sum = new double[StateoCN];
                        double sumD = 0.0;
                        
                        for (int l_d = 0; l_d < StateoD; l_d++)
                        {
                            for (int l_c = 0; l_c < StateoCN; l_c++)
                            {
                                for (int l_n = 0; l_n < StateoCN; l_n++)
                                {
                                    sum[l_n] += Math.exp(fwd[t - 1][l_d][l_c][l_n] + Math.log(D2D[l_d][k]) + Math.log(C2C[l_c][i_c]) +
                                        Math.log(Weight * N2N[l_n][i_n] + (1 - Weight) * CND[i_n][i_c + k * 2]));
                                } 
                            }
                        }

                        for (int l_n = 0; l_n < StateoCN; l_n++)
                        {
                            fwdoDCNN[t][k][i_c][i_n][l_n] = Math.log(sum[l_n]) + p;
                            sumD += sum[l_n];
                        }
                        fwd[t][k][i_c][i_n] = Math.log(sumD) + p;
                        
                        //scalingD[t] += fwd[t][k][i_c][i_n];
                    }
                }
            }



            double x = 0;

	        for (int k = 0; k < StateoD; k++)
	        {
	            for (int i_c = 0; i_c < StateoCN; i_c++)
	            {
	                for (int i_n = 0; i_n < StateoCN; i_n++)
	                {
	                	x += Math.exp(-fwd[t][0][0][0]+fwd[t][k][i_c][i_n]);
	                }
	            }
	        }


	        scalingD[t] = Math.log(x) + fwd[t][0][0][0];


            if (scalingD[t] != 0)
            {
                for (int k = 0; k < StateoD; k++)
                {
                    for (int i_c = 0; i_c < StateoCN; i_c++)
                    {
                        for (int i_n = 0; i_n < StateoCN; i_n++)
                        {
                            fwd[t][k][i_c][i_n] = fwd[t][k][i_c][i_n] - scalingD[t];
                            for (int l_n = 0; l_n < StateoCN; l_n++)
                            {
                                fwdoDCNN[t][k][i_c][i_n][l_n] = fwdoDCNN[t][k][i_c][i_n][l_n] - scalingD[t];
                            }
                        }
                    }
                }
            }
        }

        return fwd;
	}

	public void DoTwoChains(int[][] arrayC, int[][] arrayN, int iterations,
			double tolerance) throws Exception{
		if(tolerance == 0 || iterations == 0){
			throw new Exception("Iterations, tolerance should not be 0.");
		}
		
		int currentIteration = 1;
		int nc = arrayC.length;
		int nn = arrayN.length;
        boolean stop = false;
        
        double oldLikelihoodC = Double.MIN_VALUE;
        double newLikelihoodC = 0;
        double oldLikelihoodN = Double.MIN_VALUE;
        double newLikelihoodN = 0;
        
	double[][] GoC = new double[T][StateoCN];
	double[][] GoN = new double[T][StateoCN];

        do // Until convergence or max iterations
        {
            // For each sequence in the array input
            System.out.println("Current iteration: "+currentIteration);
            
            if(currentIteration == 4)
            {
            	System.out.println("Current iteration: "+currentIteration);
            }

            double[][] fwdoC = forwardC(arrayC);
            double[][] bwdoC = backwardC(arrayC);
            double[][] fwdoN = forwardN(arrayN);
            double[][] bwdoN = backwardN(arrayN);

            
			// Calculate gamma values for next computations
            for (int t = 0; t < T; t++)
            {
                double sN = 0, sC = 0, baseN = 0, baseC = 0;

                for (int k = 0; k < StateoCN; k++)
                {
                	baseC = fwdoC[t][0] + bwdoC[t][0];
                	baseN = fwdoN[t][0] + bwdoN[t][0];
                    sC += GoC[t][k] = Math.exp(fwdoC[t][k] + bwdoC[t][k] - baseC);
                    
                    sN += GoN[t][k] = Math.exp(fwdoN[t][k] + bwdoN[t][k] - baseN);
                }

                if (sC != 0 && sN != 0) // Scaling for computing
                {
                    for (int k = 0; k < StateoCN; k++)
                    {
                        GoC[t][k] /= sC;
                        GoN[t][k] /= sN;
                    }
                }
            }

            double[][][] epsilonC = new double[T][StateoCN][StateoCN];
			double[][][] epsilonN = new double[T][StateoCN][StateoCN];
			// Calculate epsilon values for next computations
            for (int t = 0; t < T - 1; t++)
            {
                double sN = 0, sC = 0;

                for (int k = 0; k < StateoCN; k++)
                {
                    for (int l = 0; l < StateoCN; l++)
                    {
                        double BTotalN = 0.0, BTotalC = 0.0;
                        for (int j = 0; j < nc; j++)
                        {
                            BTotalC = BTotalC + Math.log(C2Y[l][arrayC[j][t + 1]]);
                        }
                        for (int j = 0; j < nn; j++)
                        {
                            BTotalN = BTotalN + Math.log(N2X[l][arrayN[j][t + 1]]);
                        }
                        epsilonC[t][k][l] = fwdoC[t][k] + Math.log(C2C[k][l]) + bwdoC[t + 1][l] + BTotalC;
                        epsilonN[t][k][l] = fwdoN[t][k] + Math.log(N2N[k][l]) + bwdoN[t + 1][l] + BTotalN;
                    }
                }

                double bN = 0, bC = 0;
            	bC = epsilonC[t][0][0];
                bN = epsilonN[t][0][0];
                for (int k = 0; k < StateoCN; k++){
                    for (int l = 0; l < StateoCN; l++)
                    {
                    	epsilonC[t][k][l] = Math.exp(epsilonC[t][k][l]-bC);
                    	epsilonN[t][k][l] = Math.exp(epsilonN[t][k][l]-bN);


                    	sC += epsilonC[t][k][l];
                    	sN += epsilonN[t][k][l];
                    }
                }


                if (sC != 0 && sN != 0) // Scaling
                {
                    for (int k = 0; k < StateoCN; k++)
                        for (int l = 0; l < StateoCN; l++)
                        {
                        	epsilonC[t][k][l] /= sC;
                        	epsilonN[t][k][l] /= sN;
                        }
                }

            }

            for (int i = 0; i < scalingC.length; i++)
            {
                newLikelihoodC += scalingC[i];
                newLikelihoodN += scalingN[i];
            }

            newLikelihoodC /= nc;
            newLikelihoodN /= nn;

            System.out.println("Value C:");
            System.out.println(newLikelihoodC);
            System.out.println(C2C[0][0]);
            System.out.println(C2C[1][1]);
            System.out.println("Value N:");
            System.out.println(newLikelihoodN);
            System.out.println(N2N[0][0]);
            System.out.println(N2N[1][1]);

            // Check if the model has converged or we should stop
            if (checkConvergence(oldLikelihoodC, newLikelihoodC,
                currentIteration, iterations, tolerance) &&
                checkConvergence(oldLikelihoodN, newLikelihoodN,
                currentIteration, iterations, tolerance))
            {
                stop = true;
            }
            else
            {
                System.out.println("Continue with parameter re-estimation...");
                currentIteration++;
                oldLikelihoodC = newLikelihoodC;
                newLikelihoodC = 0.0;
                oldLikelihoodN = newLikelihoodN;
                newLikelihoodN = 0.0;

                for (int k = 0; k < StateoCN; k++)
                {
                    PioC[k] = GoC[0][k];
                    PioN[k] = GoN[0][k];

                    if (PioC[k] < 1e-10)
                    {
                        PioC[k] = 1e-10;
                    }
                    if (PioC[k] > 1 - 1e-10)
                    {
                        PioC[k] = 1 - 1e-10;
                    }

                    if (PioN[k] < 1e-10)
                    {
                        PioN[k] = 1e-10;
                    }
                    if (PioN[k] > 1 - 1e-10)
                    {
                        PioN[k] = 1 - 1e-10;
                    }
                }

                for (int i = 0; i < StateoCN; i++)
                {
                    for (int j = 0; j < StateoCN; j++)
                    {
                        double denC = 0, numC = 0, denN = 0, numN = 0;

                        for (int k = 0; k < 1; k++)
                        {
                            //int T = array[k].length;

                            for (int l = 0; l < T - 1; l++)
                            {
                                numC += epsilonC[l][i][j];
                                numN += epsilonN[l][i][j];
                            }

                            for (int l = 0; l < T - 1; l++)
                            {
                                denC += GoC[l][i];
                                denN += GoN[l][i];
                            }
                        }

                        C2C[i][j] = (denC != 0) ? numC / denC : 0.0;
                        N2N[i][j] = (denN != 0) ? numN / denN : 0.0;
                    }
                }

                for (int i = 0; i < StateoCN; i++)
                {
                    double denC = 0, denN = 0;
                    double[] numC = new double[SymboloY];
                    double[] numN = new double[SymboloX];

                    for (int k = 0; k < 1; k++)
                    {
                        for (int m = 0; m < arrayN.length; m++)
                        {
                            for (int l = 0; l < T; l++)
                            {
                                numN[arrayN[m][l]] += GoN[l][i];
                            }

                            for (int l = 0; l < T; l++)
                            {
                                denN += GoN[l][i];
                            }
                        }

                        for (int m = 0; m < arrayC.length; m++)
                        {
                            for (int l = 0; l < T; l++)
                            {
                                numC[arrayC[m][l]] += GoC[l][i];
                            }

                            for (int l = 0; l < T; l++)
                            {
                                denC += GoC[l][i];
                            }
                        }

                        for (int j = 0; j < SymboloY; j++)
                        {
                            C2Y[i][j] = (numC[j] == 0) ? 1e-10 : numC[j] / denC;
                            N2X[i][j] = (numN[j] == 0) ? 1e-10 : numN[j] / denN;
                        }
                    }
                }
            }
        } while (!stop);

	 for (int i = 0; i < 1; i++)
        {
            this.MarkC = markState(GoC);
	     this.MarkN = markState(GoN);
        }
	}

	public boolean checkConvergence(double oldLikelihood,
			double newLikelihood, int currentIteration, int maxIterations,
			double tolerance) {
		 if (tolerance > 0)
         {
             // Stopping criteria is likelihood convergence
             if (Math.abs(oldLikelihood - newLikelihood) <= tolerance)
                 return true;

             if (maxIterations > 0)
             {
                 // Maximum iterations should also be respected
                 if (currentIteration >= maxIterations)
                     return true;
             }
         }
         else
         {
             // Stopping criteria is number of iterations
             if (currentIteration == maxIterations)
                 return true;
         }

         // Check if we have reached an invalid state
         if (Double.isNaN(newLikelihood) || Double.isInfinite(newLikelihood))
         {
             return true;
         }

         return false;
	}

	public double[][] backwardN(int[][] observation) {
		double[][] bwd = new double[T][StateoCN];

        for (int i = 0; i < StateoCN; i++)
            bwd[T - 1][i] = -scalingN[T - 1];

        System.out.println("Backward...");

        for (int t = T - 2; t >= 0; t--)
        {
            for (int i = 0; i < StateoCN; i++)
            {
                double sum = 0;
                for (int j = 0; j < StateoCN; j++)
                {
                    double p = 0.0;
                    for (int k = 0; k < observation.length; k++)
                    {
                        p += Math.log(N2X[j][observation[k][t + 1]]);
                    }
                    sum += Math.exp(Math.log(N2N[i][j]) + p + bwd[t + 1][j]);
                }

                bwd[t][i] = Math.log(sum) - scalingN[t];
            }
            
        }

        return bwd;
	}

	public double[][] forwardN(int[][] observation) {
		double[][] fwd = new double[T][StateoCN];
        scalingN = new double[T];

        for (int i = 0; i < StateoCN; i++)
        {
            fwd[0][i] = Math.log(PioN[i]);
            for (int k = 0; k < observation.length; k++)
            {
                fwd[0][i] += Math.log(N2X[i][observation[k][0]]);
            }
            //scalingN[0] += fwd[0][i];
        }


        scalingN[0] = Math.log(Math.exp(-fwd[0][0]+fwd[0][1])+1) + fwd[0][0];


        if (scalingN[0] != 0) // Scaling
        {
            for (int i = 0; i < StateoCN; i++)
                fwd[0][i] = fwd[0][i] - scalingN[0];
        }

        System.out.println("Forward...");
        // 2. Induction
        for (int t = 1; t < T; t++)
        {
            for (int i = 0; i < StateoCN; i++)
            {
                double p = 0.0;
                for (int k = 0; k < observation.length; k++)
                {
                    p += Math.log(N2X[i][observation[k][t]]);
                }

                double sum = 0.0, base = 0.0;
                for (int j = 0; j < StateoCN; j++){
                	base = fwd[t - 1][0] + Math.log(N2N[0][i]);
                    sum += Math.exp((fwd[t - 1][j] + Math.log(N2N[j][i])) - base);
                }
                fwd[t][i] = Math.log(sum) + p + base;

                //scalingN[t] += fwd[t][i]; // scaling coefficient
            }

	     
            scalingN[t] = Math.log(Math.exp(-fwd[t][0]+fwd[t][1])+1) + fwd[t][0];
	     

            if (scalingN[t] != 0) // Scaling
            {
                for (int i = 0; i < StateoCN; i++)
                    fwd[t][i] = fwd[t][i] - scalingN[t];
            }
        }

        return fwd;
	}

	public double[][] backwardC(int[][] array) {
		double[][] bwd = new double[T][StateoCN];

        for (int i = 0; i < StateoCN; i++)
            bwd[T - 1][i] = -scalingC[T - 1];

        System.out.println("Backward...");
        
        //double base = 0.0;
        
        for (int t = T - 2; t >= 0; t--)
        {
            for (int i = 0; i < StateoCN; i++)
            {
                double sum = 0;
                for (int j = 0; j < StateoCN; j++)
                {
                    double p = 0;
                    for (int k = 0; k < array.length; k++)
                    {
                        p += Math.log(C2Y[j][array[k][t + 1]]);
                    }
                    sum += Math.exp((Math.log(C2C[i][j]) + p + bwd[t + 1][j]));
                }

                bwd[t][i] += Math.log(sum) - scalingC[t];
                
            }
        }

        return bwd;
	}

	public double[][] forwardC(int[][] array) {
		double[][] fwd = new double[T][StateoCN];
        scalingC = new double[T];

        for (int i = 0; i < StateoCN; i++)
        {

            fwd[0][i] = Math.log(PioC[i]);
            for (int k = 0; k < array.length; k++)
            {
                fwd[0][i] += Math.log(C2Y[i][array[k][0]]);
            }
            //scalingC[0] += fwd[0][i];
        }


        scalingC[0] = Math.log(Math.exp(-fwd[0][0]+fwd[0][1])+1) + fwd[0][0];


        if (scalingC[0] != 0) // Scaling
        {
            for (int i = 0; i < StateoCN; i++)
                fwd[0][i] = fwd[0][i] - scalingC[0];
        }

        System.out.println("Forward...");
        // 2. Induction
        for (int t = 1; t < T; t++)
        {
            for (int i = 0; i < StateoCN; i++)
            {
                double p = 0.0;
                for (int k = 0; k < array.length; k++)
                {
                    p += Math.log(C2Y[i][array[k][t]]);
                }

                double sum = 0.0, base = 0.0;
                for (int j = 0; j < StateoCN; j++){
                	base = fwd[t - 1][0] + Math.log(C2C[0][i]);
                    sum += Math.exp((fwd[t - 1][j] + Math.log(C2C[j][i])) - base);
                }
                fwd[t][i] = Math.log(sum) + p + base;

                //scalingC[t] += fwd[t][i]; // scaling coefficient
                
            }


            scalingC[t] = Math.log(Math.exp(-fwd[t][0]+fwd[t][1])+1) + fwd[t][0];


            if (scalingC[t] != 0) // Scaling
            {
                for (int i = 0; i < StateoCN; i++)
                    fwd[t][i] = fwd[t][i] - scalingC[t];
            }
        }

        return fwd;
	}

}
