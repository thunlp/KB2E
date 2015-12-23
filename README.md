New: Add PTransE (EMNLP 2015) code!



==RESULT==

We list the result of the code in date set FB15k and WN18.

FB15k

| Model      |     MeanRank(Raw) |   MeanRank(Filter)   |	Hit@10(Raw)	| Hit@10(Filter)|
| :-------- | --------:| :------: | :------: |:------: |
| TransE(paper)|    243 | 125 |  34.9 | 47.1|
| TransH(paper)        |   212 |  87 |  45.7 | 64.4|
| TransR(n=50)        |    198| 77 |  48.2 | 68.7 |
| TransE(Our, n=50)   | 210|	82  |	41.9|  61.3 |
| TransE(Our, n=100)  |    205 |  63 |  47.9 | 70.2 |
|PTransE (ADD, 2-step) |    200 | 54 | 51.8 | 83.4|
|PTransE (MUL, 2-step) |    216 |  67 | 47.4 | 77.7 |
|PTransE (RNN, 2-step) | 242 | 92 | 50.6 | 82.2 |
|PTransE (ADD, 3-step) |207 | 58 | 51.4 | 84.6 |

WN18

| Model      |     MeanRank(Raw) |   MeanRank(Filter)   |	Hit@10(Raw)	| Hit@10(Filter)|
| :-------- | --------:| :------: | :------: |:------: |
| TransE(paper)|    263 |    251 | 75.4 | 89.2|
| TransH(paper)        |    318 |    303 | 75.4 | 86.7|
| TransR        |    238 | 225 | 79.8  |92.0|
| TransE(Our)   | 251	|239|78.9|		89.8|


===== DATA =====

I provide FB15k and WN18  datasets used for the task link prediction with the input format of my code in data.zip.

The original data use in the experiment can download in:

FB15k, WN18 are published by the author of the paper "Translating Embeddings for Modeling Multi-relational Data (2013)." [[Download]](https://everest.hds.utc.fr/doku.php?id=en:transe)

FB13, WN11 are published by the author of the paper "Reasoning With Neural Tensor Networks for Knowledge Base Completion". [[Download]](http://cs.stanford.edu/~danqi/data/nips13-dataset.tar.bz2)

New York Times Corpus:  The data used in relation extraction from text which is publish by the paper " Modeling relations and their mentions without labeled text". If you want the data, you should buy from LDC (https://catalog.ldc.upenn.edu/LDC2008T19) first.

FB40k [[Download]](http://pan.baidu.com/s/1c0xrtVa)



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

	call the program Test_TransR with method as parameter in folder TransR/

CTransR:

	call the program Test_CTransR with method as parameter in folder CTransR/

It will evaluate on test.txt and report mean rank and Hits@10




==CITE==

If you use the code, you should cite the following paper:

Yankai Lin, Zhiyuan Liu, Maosong Sun, Yang Liu, Xuan Zhu. Learning Entity and Relation Embeddings for Knowledge Graph Completion. The 29th AAAI Conference on Artificial Intelligence (AAAI'15).[[pdf]](http://nlp.csai.tsinghua.edu.cn/~lzy/publications/aaai2015_transr.pdf)
