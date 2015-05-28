#include<iostream>
#include<cstring>
#include<cstdio>
#include<map>
#include<vector>
#include<string>
#include<ctime>
#include<algorithm>
#include<cmath>
#include<cstdlib>
#include<pthread.h>
using namespace std;

bool debug=false;
bool L1_flag=1;

string version;
string trainortest = "test";

map<string,int> relation2id,entity2id;
map<int,string> id2entity,id2relation;
map<string,string> mid2name,mid2type;
map<int,map<int,int> > entity2num;
map<int,int> e2num;
map<pair<string,string>,map<string,double> > rel_left,rel_right;

int relation_num,entity_num,fb_relation_num;
map<int,int> cluster2fb;
map<int,vector<int> > fb2cluster;
int n = 100;
int m = 100;

double sigmod(double x)
{
    return 1.0/(1+exp(-x));
}

double vec_len(vector<double> a)
{
	double res=0;
	for (int i=0; i<a.size(); i++)
		res+=a[i]*a[i];
	return sqrt(res);
}

void vec_output(vector<double> a)
{
	for (int i=0; i<a.size(); i++)
	{
		cout<<a[i]<<"\t";
		if (i%9==4)
			cout<<endl;
	}
	cout<<"-------------------------"<<endl;
}

double sqr(double x)
{
    return x*x;
}

char buf[100000],buf1[100000];

int my_cmp(pair<double,int> a,pair<double,int> b)
{
    return a.first>b.first;
}

double cmp(pair<int,double> a, pair<int,double> b)
{
	return a.second<b.second;
}

class Test{
    vector<vector<double> > relation_vec,entity_vec,fb_relation_vec;
     vector<vector<vector<double> > > A, fb_A;


    vector<int> h,l,r;
    vector<int> fb_h,fb_l,fb_r;
    map<pair<int,int>, map<int,int> > ok;
    double res ;
public:
    void add(int x,int y,int z, bool flag)
    {
    	if (flag)
    	{
        	fb_h.push_back(x);
        	fb_r.push_back(z);
        	fb_l.push_back(y);
        }
        ok[make_pair(x,z)][y]=1;
    }

    int rand_max(int x)
    {
        int res = (rand()*rand())%x;
        if (res<0)
            res+=x;
        return res;
    }
    double len;
    double calc_sum(int e1,int e2,int rel)
    {
        double cluster_rate = 1;
    	int fb_rel = cluster2fb[rel];
    	vector<double> e1_vec;
        e1_vec.resize(m);
        vector<double> e2_vec;
        e2_vec.resize(m);
        for (int ii=0; ii<m; ii++)
        {
            for (int jj=0; jj<n; jj++)
            {
                e1_vec[ii]+=A[fb_rel][jj][ii]*entity_vec[e1][jj];
                e2_vec[ii]+=A[fb_rel][jj][ii]*entity_vec[e2][jj];
            }
        }
        double sum=0;
        if (L1_flag)
        	for (int ii=0; ii<m; ii++)
            	sum+=-fabs(e2_vec[ii]-e1_vec[ii]-cluster_rate*relation_vec[rel][ii]-(1-cluster_rate)*fb_relation_vec[fb_rel][ii]);
        else
        	for (int ii=0; ii<m; ii++)
            	sum+=-sqr(e2_vec[ii]-e1_vec[ii]-cluster_rate*relation_vec[rel][ii]-(1-cluster_rate)*fb_relation_vec[fb_rel][ii]);
        return sum;
    }
    void run()
    {

        FILE* f1 = fopen(("fb_relation2vec"+version).c_str(),"r");
        FILE* f2 = fopen(("relation2vec"+version).c_str(),"r");
        FILE* f3 = fopen(("entity2vec"+version).c_str(),"r");
        cout<<relation_num<<' '<<entity_num<<' '<<fb_relation_num<<endl;
        fb_relation_vec.resize(fb_relation_num);
        for (int i=0; i<relation_num;i++)
        {
            fb_relation_vec[i].resize(n);
            for (int ii=0; ii<n; ii++)
                fscanf(f1,"%lf",&fb_relation_vec[i][ii]);
        }
        relation_vec.resize(relation_num);
        for (int i=0; i<relation_num;i++)
        {
            relation_vec[i].resize(n);
            for (int ii=0; ii<n; ii++)
                fscanf(f2,"%lf",&relation_vec[i][ii]);
        }

        entity_vec.resize(entity_num);
        cout<<entity_num<<endl;
        for (int i=0; i<entity_num;i++)
        {
            entity_vec[i].resize(n);
            for (int ii=0; ii<n; ii++)
                fscanf(f3,"%lf",&entity_vec[i][ii]);
        }

        FILE* f4 = fopen(("A.txt"+version).c_str(),"r");

        A.resize(fb_relation_num);
		for (int i=0; i<fb_relation_num; i++)
		{
		    A[i].resize(n);
		    for (int jj=0; jj<n; jj++)
		    {
		        A[i][jj].resize(m);
		        for (int ii=0; ii<m; ii++)
                    fscanf(f4,"%lf",&A[i][jj][ii]);
		    }
		}
		fclose(f4);
        cout<<"Here"<<endl;
		double lsum=0 ,lsum_filter= 0;
		double rsum = 0,rsum_filter=0;
		double lp_n=0,lp_n_filter=0;
		double rp_n=0,rp_n_filter=0;
		map<int,double> lsum_r,lsum_filter_r;
		map<int,double> rsum_r,rsum_filter_r;
		map<int,double> lp_n_r,lp_n_filter_r;
		map<int,double> rp_n_r,rp_n_filter_r;
		map<int,int> rel_num;
        //cout<<"here"<<endl;


        FILE* f_test_cluster = fopen("../cluster/output/test_cluster.txt","r");
        for (int testid = 0; testid<fb_l.size(); testid+=1)
		{
		    vector<int> h_rel,l_rel;
		    for (int i=0; i<entity_num; i++)
            {
                int x;
                fscanf(f_test_cluster,"%d",&x);
                h_rel.push_back(x);
            }
            for (int i=0; i<entity_num; i++)
            {
                int x;
                fscanf(f_test_cluster,"%d",&x);
                l_rel.push_back(x);
            }
			int h = fb_h[testid];
			int l = fb_l[testid];
			int rel = fb_r[testid];



            if (h_rel[l]!=l_rel[h])
                cout<<"wrong"<<h_rel[l]<<' '<<l_rel[h]<<' '<<h<<' '<<l<<endl;

			double tmp = 0;
			rel_num[rel]+=1;
			vector<pair<int,double> > a;
			for (int i=0; i<entity_num; i++)
			{
				double sum  = calc_sum(i,l,h_rel[i]);
			//	for (int j=0; j<fb2cluster[rel].size(); j++)
              //      sum+= calc_sum(i,l,fb2cluster[rel][j]);
				a.push_back(make_pair(i,sum));
			}
			sort(a.begin(),a.end(),cmp);
			double ttt=0;
			int filter = 0;
			for (int i=a.size()-1; i>=0; i--)
			{
			    if (ok[make_pair(a[i].first,rel)].count(l)==0)
			    	filter+=1;
				if (a[i].first ==h)
				{
					lsum+=a.size()-i;
					lsum_filter+=filter+1;
					lsum_r[rel]+=a.size()-i;
					lsum_filter_r[rel]+=filter+1;
					if (a.size()-i<=10)
					{
						lp_n+=1;
						lp_n_r[rel]+=1;
					}
					if (filter<10)
					{
						lp_n_filter+=1;
						lp_n_filter_r[rel]+=1;
					}
					break;
				}
			}
			a.clear();
			for (int i=0; i<entity_num; i++)
			{
				double sum =  calc_sum(h,i,l_rel[i]);
				//for (int j=0; j<fb2cluster[rel].size(); j++)
                  //  sum+= calc_sum(h,i,fb2cluster[rel][j]);
				a.push_back(make_pair(i,sum));
			}
			sort(a.begin(),a.end(),cmp);
			ttt=0;
			filter=0;
			for (int i=a.size()-1; i>=0; i--)
			{
			    if (ok[make_pair(h,rel)].count(a[i].first)==0)
			    	filter+=1;
				if (a[i].first==l)
				{
					rsum+=a.size()-i;
					rsum_filter+=filter+1;
					rsum_r[rel]+=a.size()-i;
					rsum_filter_r[rel]+=filter+1;
					if (a.size()-i<=10)
					{
						rp_n+=1;
						rp_n_r[rel]+=1;
					}
					if (filter<10)
					{
						rp_n_filter+=1;
						rp_n_filter_r[rel]+=1;
					}
					break;
				}
			}
			//if (testid%100==0)
			//cout<<testid<<":"<<"\t"<<lsum/(testid+1)<<' '<<lp_n/(testid+1)<<' '<<rsum/(testid+1)<<' '<<rp_n/(testid+1)<<"\t"<<lsum_filter/(testid+1)<<' '<<lp_n_filter/(testid+1)<<' '<<rsum_filter/(testid+1)<<' '<<rp_n_filter/(testid+1)<<endl;
		}
		cout<<"left:"<<lsum/fb_l.size()<<'\t'<<lp_n/fb_l.size()<<"\t"<<lsum_filter/fb_l.size()<<'\t'<<lp_n_filter/fb_l.size()<<endl;
		cout<<"right:"<<rsum/fb_r.size()<<'\t'<<rp_n/fb_r.size()<<'\t'<<rsum_filter/fb_r.size()<<'\t'<<rp_n_filter/fb_r.size()<<endl;
		for (int rel=0; rel<fb_relation_num; rel++)
		{
			int num = rel_num[rel];
			cout<<"rel:"<<id2relation[rel]<<' '<<num<<endl;
			cout<<"left:"<<lsum_r[rel]/num<<'\t'<<lp_n_r[rel]/num<<"\t"<<lsum_filter_r[rel]/num<<'\t'<<lp_n_filter_r[rel]/num<<endl;
			cout<<"right:"<<rsum_r[rel]/num<<'\t'<<rp_n_r[rel]/num<<'\t'<<rsum_filter_r[rel]/num<<'\t'<<rp_n_filter_r[rel]/num<<endl;
		}
    }

};
Test test;

void prepare()
{
    FILE* f1 = fopen("../data/entity2id.txt","r");
	FILE* f2 = fopen("../data/relation2id.txt","r");
	int x;
	while (fscanf(f1,"%s%d",buf,&x)==2)
	{
		string st=buf;
		entity2id[st]=x;
		id2entity[x]=st;
		mid2type[st]="None";
		entity_num++;
	}
	while (fscanf(f2,"%s%d",buf,&x)==2)
	{
		string st=buf;
		relation2id[st]=x;
		id2relation[x]=st;
		fb_relation_num++;
	}


	FILE* f_kb_train = fopen("../data/train.txt","r");
	FILE* f_train_cluster = fopen("../cluster/output/train_cluster.txt","r");
	while (fscanf(f_kb_train,"%s",buf)==1)
    {
        string s1=buf;
        fscanf(f_kb_train,"%s",buf);
        string s2=buf;
        fscanf(f_kb_train,"%s",buf);
        string s3=buf;
        int rel_fb = relation2id[s3];
        int rel;
        fscanf(f_train_cluster,"%d",&rel);
        if (cluster2fb.count(rel)==0)
        {
            relation_num+=1;
            cluster2fb[rel]=rel_fb;
            fb2cluster[rel_fb].push_back(rel);
        }
        else
        if (cluster2fb[rel]!=rel_fb)
        {
            cout<<"not map:"<<entity2id[s1]<<" "<<entity2id[s2]<<" "<<rel<<' '<<rel_fb<<endl;
            exit(-1);
        }
        test.add(entity2id[s1],entity2id[s2],relation2id[s3],false);
    }

    fclose(f_train_cluster);
    fclose(f_kb_train);

    FILE* f_kb_test = fopen("../data/test.txt","r");
	while (fscanf(f_kb_test,"%s",buf)==1)
    {
        string s1=buf;
        fscanf(f_kb_test,"%s",buf);
        string s2=buf;
        fscanf(f_kb_test,"%s",buf);
        string s3=buf;
        test.add(entity2id[s1],entity2id[s2],relation2id[s3],true);
    }
    fclose(f_kb_test);

    FILE* f_kb_valid = fopen("../data/valid.txt","r");
	while (fscanf(f_kb_valid,"%s",buf)==1)
    {
        string s1=buf;
        fscanf(f_kb_valid,"%s",buf);
        string s2=buf;
        fscanf(f_kb_valid,"%s",buf);
        string s3=buf;
        test.add(entity2id[s1],entity2id[s2],relation2id[s3],false);
    }
    fclose(f_kb_valid);
}


int main(int argc,char**argv)
{
    if (argc<2)
        return 0;
    else
    {
        version = argv[1];
        if (argc==3)
            trainortest = argv[2];
        prepare();
        test.run();
    }
}

