#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;


map<string,int> relation2id,entity2id;
map<int,string> id2entity,id2relation;
int entity_num=0,relation_num=0;
char buf[100000],buf1[100000];
vector<vector<double> > entity_vec;
vector<double> a;

double sqr(double x)
{
    return x*x;
}


void vec_output(vector<double> &a)
{
    for (int i=0; i<a.size(); i++)
        cout<<a[i]<<'\t';
    cout<<endl;
}

int main(int argc, char** argv)
{
    int rel_now=-1;
    if (argc ==2)
        rel_now =  atoi(argv[1]);
    else
        return 0;
    cout<<rel_now<<endl;


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
    FILE* f_kb = fopen("../data/train.txt","r");
    vector<int> train_e1,train_e2,train_rel,train_index;
    int index=-1;
	while (fscanf(f_kb,"%s",buf)==1)
    {
        string s1=buf;
        fscanf(f_kb,"%s",buf);
        string s2=buf;
        fscanf(f_kb,"%s",buf);
        string s3=buf;
        int e1 = entity2id[s1];
        int e2 = entity2id[s2];
        int rel = relation2id[s3];
        index+=1;
        if (rel!=rel_now)
            continue;
        train_index.push_back(index);
        train_e1.push_back(e1);
        train_e2.push_back(e2);
        train_rel.push_back(rel);
    }
    FILE* f_ent_vec = fopen("../TransE/entity2vec.txt2","r");

   // FILE* f_ent_vec = fopen("../data/entity_vec_word2vec.txt","r");
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
    double sum = 0;
    int n=train_e1.size();
    FILE* g_sim = fopen("similarity.txt","w");
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
        if (i!=j)
        {
            double tmp=0;
            for (int ii=0; ii<vec_size; ii++)
                tmp-=sqr(entity_vec[train_e1[i]][ii]-entity_vec[train_e2[i]][ii]-
                          entity_vec[train_e1[j]][ii]+entity_vec[train_e2[j]][ii]);
            a.push_back(tmp);
            if (tmp>1e6)
            {
                vec_output(entity_vec[train_e1[i]]);
                vec_output(entity_vec[train_e2[i]]);
                vec_output(entity_vec[train_e1[j]]);
                vec_output(entity_vec[train_e2[j]]);
                exit(0);
            }
            fprintf(g_sim,"%d\t%d\t%lf\n",i+1,j+1,tmp);
        }
    fclose(g_sim);
    if (n>1)
    {
        sort(a.begin(),a.end());
        cout<<a[0]<<' '<<a[a.size()-1]<<' '<<a[0]-(a[a.size()-1]-a[0])<<endl;
        sum = a[0]-10*(a[a.size()-1]-a[0]);
    }
    FILE* g_pre = fopen("pre.txt","w");
    for (int i=0; i<n; i++)
        fprintf(g_pre,"%lf\n",sum);
    fclose(g_pre);
    FILE* g_env = fopen(("data/id2index"+string(argv[1])+".txt").c_str(),"w");
    for (int i=0; i<n; i++)
    {
        fprintf(g_env,"%d\n",train_index[i]);
    }
    fclose(g_env);
    return 0;
}
