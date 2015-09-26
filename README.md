# Collaborative-Filtering
The project uses Movielens dataset (the 100K dataset) which is available open source.
It computes the following four algorithms (on the raw ratings data, without baseline correction) and calculate how the NMAE (Normalized Mean Absolute Error) varies.
1. Non-Negative Matrix Factorization (NNMF)
2. Incremented Rank Power Factorization (IRPF)
3. Fixed Point Continuation (FPC)
4. Singular Value Thresholding (SVT)
For computing NMAE in case of IRPF, a globally available function licensed under the CC Attribution-Noncommercial-Share Alike 
% 3.0 is used.
The evaluation metric uses 5-fold cross validation for reporting the results. 
