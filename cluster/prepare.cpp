#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>

using namespace std;


map<string,int> relation2id,entity2id;
map<int,string> id2entity,id2relation;
int entity_num=0,relation_num=0;
char buf[100000],buf1[100000];
vector<vector<double> > entity_vec;
vector<vector<int> > rel_cluster_set;


double sqr(double x)
{
    return x*x;
}

vector<vector<double> > cluster_center;
vector<int> cluster_size;

int main()
{
    FILE* f1 = fopen("../data/entity2id.txt","r");
	FILE* f2 = fopen("../data/relation2id.txt","r");
	int x;
	while (fscanf(f1,"%s%d",buf,&x)==2)
	{
		string st=buf;
		entity2id[st]=x;
		id2entity[x]=st;
		entity_num++;
	}
	while (fscanf(f2,"%s%d",buf,&x)==2)
	{
		string st=buf;
		relation2id[st]=x;
		id2relation[x]=st;
		relation_num++;
	}
	fclose(f1);
	fclose(f2);
    FILE* f_ent_vec = fopen("../TransE/entity2vec.txt2","r");
    int vec_size=50;
    entity_vec.resize(entity_num);
    for (int i=0; i<entity_num;i++)
    {
        entity_vec[i].resize(vec_size);
        for (int ii=0; ii<vec_size; ii++)
        {
            fscanf(f_ent_vec,"%lf",&entity_vec[i][ii]);
        }
    }
    fclose(f_ent_vec);
    FILE* f_kb_train = fopen("../data/train.txt","r");
    vector<int> train_e1,train_e2,train_rel,train_index,train_cluster;
    int train_num=0;
	while (fscanf(f_kb_train,"%s",buf)==1)
    {
        string s1=buf;
        fscanf(f_kb_train,"%s",buf);
        string s2=buf;
        fscanf(f_kb_train,"%s",buf);
        string s3=buf;
        int e1 = entity2id[s1];
        int e2 = entity2id[s2];
        int rel = relation2id[s3];
        train_num+=1;
        train_e1.push_back(e1);
        train_e2.push_back(e2);
        train_rel.push_back(rel);
    }
    train_cluster.resize(train_num);
    fclose(f_kb_train);
    int cluster_num=0;
    map<int,int> cluster2id;
    rel_cluster_set.resize(relation_num);
    for (int rel =0; rel<relation_num; rel++)
    {
        map<int,int> ok;
        ok.clear();
        string s;
        stringstream ss;
        ss<<rel;
        ss>>s;
        cout<<s<<endl;
        FILE* f_index = fopen(("data/id2index"+s+".txt").c_str(),"r");
        FILE* f_cluster = fopen(("data/cluster"+s+".txt").c_str(),"r");
        int x,y;
        while (fscanf(f_index,"%d",&x)==1)
        {
            if (f_cluster!=NULL)
                fscanf(f_cluster,"%d",&y);
            else
                y=0;
            if (ok.count(y)==0)
            {
                ok[y]=cluster_num;
                cluster2id[cluster_num]=x;
                rel_cluster_set[train_rel[x]].push_back(cluster_num);

                cluster_size.push_back(0);
                vector<double> tmp;
                tmp.resize(vec_size);
                cluster_center.push_back(tmp);
                cluster_num+=1;
            }
            for (int ii=0; ii<vec_size; ii++)
                cluster_center[ok[y]][ii]+=entity_vec[train_e1[x]][ii]-entity_vec[train_e2[x]][ii];
            train_cluster[x]=ok[y];
            cluster_size[ok[y]]+=1;
        }
        fclose(f_index);
        if (f_cluster!=NULL)
            fclose(f_cluster);
    }
    FILE* f_train_cluster = fopen("output/train_cluster.txt","w");
    for (int i=0; i<train_num; i++)
        fprintf(f_train_cluster,"%d\n",train_cluster[i]);
    fclose(f_train_cluster);
    cout<<cluster_num<<endl;
    FILE* f_test_cluster = fopen("output/test_cluster.txt","w");
    FILE* f_kb_test = fopen("../data/test.txt","r");
   // for (int i=0; i<cluster_num; i++)
   // {
  //      cout<<cluster_center[i][0]/cluster_size[i]<<endl;
   // }
    int step=0;
	while (fscanf(f_kb_test,"%s",buf)==1)
    {
        cout<<step++<<endl;
        string s1=buf;
        fscanf(f_kb_test,"%s",buf);
        string s2=buf;
        fscanf(f_kb_test,"%s",buf);
        string s3=buf;
        int e1 = entity2id[s1];
        int e2 = entity2id[s2];
        int rel = relation2id[s3];
       // cout<<"here"<<endl;
        for (int h=0; h<entity_num; h++)
        {
            double dis = 1e9;
            int k=-1;
            for (int i=0; i<rel_cluster_set[rel].size(); i++)
            {
                double tmp = 0;
                int x = rel_cluster_set[rel][i];
               // cout<<x<<' '<<rel_cluster_set[rel][i]<<' '<<h<<' '<<train_e1[x]<<' '<<train_e2[x]<<endl;
                for (int ii=0; ii<50; ii++)
                    tmp+=sqr(entity_vec[h][ii]-entity_vec[e2][ii]-cluster_center[x][ii]/cluster_size[x]);

               // cout<<x<<endl;
                if (tmp<dis)
                {
                    dis = tmp;
                    k = rel_cluster_set[rel][i];
                }
            }
            fprintf(f_test_cluster,"%d\t",k);
        }

       // cout<<"here"<<endl;
        fprintf(f_test_cluster,"\n");
        for (int t=0; t<entity_num; t++)
        {
            double dis = 1e9;
            int k=-1;
            for (int i=0; i<rel_cluster_set[rel].size(); i++)
            {
                double tmp = 0;
                int x = rel_cluster_set[rel][i];
                for (int ii=0; ii<50; ii++)
                    tmp+=sqr(entity_vec[e1][ii]-entity_vec[t][ii]-cluster_center[x][ii]/cluster_size[x]);
               // cout<<"tmp="<<tmp<<' '<<rel_cluster_set[rel][i]<<endl;
                if (tmp<dis)
                {
                    dis = tmp;
                    k = rel_cluster_set[rel][i];
                }
            }
            fprintf(f_test_cluster,"%d\t",k);
        }
        fprintf(f_test_cluster,"\n");
    }
    fclose(f_kb_test);
    return 0;
}
