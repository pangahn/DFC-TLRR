# DFC-TLRR

[![HitCount](http://hits.dwyl.io/ploxolg/DFC-TLRR.svg)](http://hits.dwyl.io/ploxolg/DFC-TLRR)

The MATLAB implementation for distributed low-rank subspace clustering (DFC-TLRR) method for tensor data, proposed in Distributed Low-rank Tensor Subspace Clustering.

Accepted by Journal of South China University of Technology (Natural Science Edition) 2019.02.

You should first add the [tensor toolbox](https://github.com/andrewssobral/tensor_toolbox) into you matlab path.
    
> Abstract: Low-Rank Representation (LRR), a state-of-the-art method for subspace clustering problem, is provably and empirically accurate on small problems but cannot handle large-scale problems effectively. Distributed low-rank subspace clustering (DFC-LRR) introduces a novel divide-and-conquer algorithm for large-scale subspace clustering problems, but it cannot handle the multi-way data directly. To solve this issue, we propose a distributed low-rank subspace clustering algorithm (DFC-TLRR) for tensor data. The proposed method considers multi-way data as tensor and extends LRR subspace clustering method to multi-way data by combining t-product operator. And it adopts the distributed parallel computing to obtain the low-rank coefficient tensor, then get the sparse similarity matrix by sparing every lateral slices of the coefficient tensor. Experimental results on the Extended Yale B,COIL20 and UCSD datasets show that the proposed algorithm outperforms DFC-LRR in clustering accuracy, and distributed computing can reduce the running time obviously.
Key words: low-rank representation; subspace clustering; distributed computing; tensor
