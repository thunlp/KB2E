===== DATA =====

The data use in the experiment can download in:

FB15k, WN18 are published by the author of the paper "Translating Embeddings for Modeling Multi-relational Data (2013)."

FB13, WN11 are published by the author of the paper "Reasoning With Neural Tensor Networks for Knowledge Base Completion".



Datasets are needed in the folder data/ in the following format

Dataset contains six files:



+ train.txt: training file, format (e1, e2, rel).

+ valid.txt: validation file, same format as train.txt

+ test.txt: test file, same format as train.txt.

+ entity2id.txt: all entities and corresponding ids, one per line.

+ relation2id.txt: all relations and corresponding ids, one per line.



Currently we cannot upload data due to huge size. We will release data with codes together once the paper is published.



===== CODE =====

In the folder TransE/, TransR/, CTransR/:



===== COMPILE =====

Just type make in the folder ./



== TRAINING ==

For training, You need follow the step below:





TransE:

	call the program Train_TransE in folder TransE/
	
TransH:
	call the program Train_TransH in folder TransH/

TransR:

	1:	Train the unif method of TransE as initialization.

	2:  call the program Train_TransR in folder TransR/

CTransR:

	1:	Train the unif method of TransR as initialization.

	2:  run the bash run.sh with relation number in folder cluster/ to cluster the triples in the trainning data.

		i.e.

			bash run.sh 10

	3:  call the program Train_cTransR in folder CTransR/

You can also change the parameters when running Train_TransE, Train_TransR, Train_CTransR.

-size : the embedding size k, d

-rate : learing rate

-method: 0 - unif, 1 - bern



== TESTING ==

For testing, You need follow the step below:


TransR:

	call the program Train_TransR with method as parameter in folder TransR/

CTransR:

	call the program Train_CTransR with method as parameter in folder CTransR/

It will evaluate on test.txt and report mean rank and Hits@10

==CITE==

If you use the code, you should cite the following paper:

Lin, Y., Liu, Z., Sun, M., Liu, Y., & Zhu, X. (2015). Learning Entity and Relation Embeddings for Knowledge Graph Completion.

