## The project will no longer be maintained and users are recommended to access and use the new package https://github.com/thunlp/OpenKE.

Evaluation Results
==========

We list the result of various methods implemented by ourselves in dateset FB15k and WN18.

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


Data
==========

We provide FB15k and WN18 datasets used for the task link prediction in data.zip, using the input format required by our codes. The original data can be downloaded from:

FB15k, WN18 are published by "Translating Embeddings for Modeling Multi-relational Data (2013)." [[Download]](https://everest.hds.utc.fr/doku.php?id=en:transe)

FB13, WN11 are published by "Reasoning With Neural Tensor Networks for Knowledge Base Completion". [[Download]](http://cs.stanford.edu/~danqi/data/nips13-dataset.tar.bz2)

New York Times Corpus: The data used in relation extraction from text is publish by "Modeling relations and their mentions without labeled text". The data should be obtained from LDC (https://catalog.ldc.upenn.edu/LDC2008T19) first.

FB40k [[Download]](http://pan.baidu.com/s/1c0xrtVa)

Datasets are required in the folder data/ in the following format, containing six files:

+ train.txt: training file, format (e1, e2, rel).

+ valid.txt: validation file, same format as train.txt

+ test.txt: test file, same format as train.txt.

+ entity2id.txt: all entities and corresponding ids, one per line.

+ relation2id.txt: all relations and corresponding ids, one per line.

+ e1_e2.txt: the top-500 entity pairs which are calculated by TransE.  [[Download]](https://pan.baidu.com/s/1c2iLtmg)

Code
==========

The codes are in the folder TransE/, TransR/, CTransR/.

Compile
==========

Just type make in the folder ./

Training
==========

For training, you need to follow the steps below:

TransE: call the program Train_TransE in folder TransE/
	
TransH: call the program Train_TransH in folder TransH/

TransR:

+ Train the unif method of TransE as initialization.

+ Call the program Train_TransR in folder TransR/

CTransR:

+ Train the unif method of TransR as initialization.

+ Run the bash run.sh with relation number in folder cluster/ to cluster triples in the trainning data. e.g., bash run.sh 10

+ Call the program Train_cTransR in folder CTransR/

You can also change the parameters when running Train_TransE, Train_TransR, Train_CTransR.

-size : the embedding size k, d

-rate : learing rate

-method: 0 - unif, 1 - bern

Testing
==========

For testing, you need to follow the steps below:

TransR: Call the program Test_TransR with method as parameter in folder TransR/

CTransR: Call the program Test_CTransR with method as parameter in folder CTransR/

It will evaluate on test.txt and report mean rank and Hits@10.

Cite
==========

If you use the code, please kindly cite the following paper:

Yankai Lin, Zhiyuan Liu, Maosong Sun, Yang Liu, Xuan Zhu. Learning Entity and Relation Embeddings for Knowledge Graph Completion. The 29th AAAI Conference on Artificial Intelligence (AAAI'15).[[pdf]](http://nlp.csai.tsinghua.edu.cn/~lzy/publications/aaai2015_transr.pdf)
