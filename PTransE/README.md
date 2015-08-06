===== DATA =====

I provide FB15k  datasets used for the task Knowledge Base Completion with the input format in  [[Download]](http://pan.baidu.com/s/1qWwuP1U)

Dataset contains six files:



+ train.txt: training file, format (e1, e2, rel).

+ valid.txt: validation file, same format as train.txt

+ test.txt: test file, same format as train.txt.

+ entity2id.txt: all entities and corresponding ids, one per line.

+ relation2id.txt: all relations and corresponding ids, one per line.

+ e1_e2.txt: all top-500 entity pairs mentioned in the task entity prediction.



===== CODE =====

In the folder PTransE_add/, PTransE_mul/, PTransE_RNN/:



===== COMPILE =====

Just type make in the folder ./



== TRAINING ==

For training, You need follow the step below:

1. Data preprocessing: python PCRA.py

2. call Train_TransE_path in the corresponding directory.
	./Train_TransE_path 1


== TESTING ==

For testing, You need follow the step below:

call Test_TransE_path in the corresponding directory.
	./Test_TransE_path 1

It will evaluate on test.txt and report mean rank and Hits@1, the format is 
tier: head_mean_rank(Raw) head_hit@1(Raw), head_mean_rank(Filter) head_hit@1(Filter), tail_mean_rank(Raw) tail_hit@1(Raw), tail_mean_rank(Filter) tail_hit@1(Filter)
relation_mean_rank(Raw) relation_hit@1(Raw), relation_mean_rank(Filter) relation_hit@1(Filter)


