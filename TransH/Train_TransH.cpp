
#include<iostream>
#include<cstring>
#include<cstdio>
#include<map>
#include<vector>
#include<string>
#include<ctime>
#include<cmath>
#include<cstdlib>
using namespace std;


#define pi 3.1415926535897932384626433832795

//normal distribution
double rand(double min, double max)
{
    return min+(max-min)*rand()/(RAND_MAX+1.0);
}
double normal(double x, double miu,double sigma)
{
    return 1.0/sqrt(2*pi)/sigma*exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}
double randn(double miu,double sigma, double min ,double max)
{
    double x,y,dScope;
    do{
        x=rand(min,max);
        y=normal(x,miu,sigma);
        dScope=rand(0.0,normal(miu,miu,sigma));
    }while(dScope>y);
    return x;
}

double sqr(double x)
{
    return x*x;
}

double vec_len(vector<double> &a)
{
	double res=0;
	for (int i=0; i<a.size(); i++)
		res+=a[i]*a[i];
	return sqrt(res);
}

string version;
char buf[100000],buf1[100000];
int relation_num,entity_num,feature_num,nyt_relation_num;
map<string,int> relation2id,entity2id;
map<int,string> id2entity,id2relation;



map<int,map<int,vector<int> > > left_entity,right_entity;
map<int,double> left_mean,right_mean,left_var,right_var;
map<int,map<int,int > > left_candidate_ok,right_candidate_ok;
map<int,vector<int> > left_candidate,right_candidate;
map<int,int> entity2num;

class Train{

public:
	map<pair<int,int>, map<int,int> > ok;
    void add(int x,int y,int z)
    {
        fb_h.push_back(x);
        fb_r.push_back(z);
        fb_l.push_back(y);
        ok[make_pair(x,z)][y]=1;
    }
    void run(int n_in,double rate_in,double margin_in,int method_in)
    {
        n = n_in;
        rate = rate_in;
        margin = margin_in;
        method = method_in;
		A.resize(relation_num);
		for (int i=0; i<relation_num; i++)
		{
		    A[i].resize(n);
		    for (int j=0; j<n; j++)
		    	A[i][j] = randn(0,1.0/n,-1,1);
		    norm2one(A[i]);
		}
        relation_vec.resize(relation_num);
		for (int i=0; i<relation_vec.size(); i++)
			relation_vec[i].resize(n);
        entity_vec.resize(entity_num);
		for (int i=0; i<entity_vec.size(); i++)
			entity_vec[i].resize(n);
        for (int i=0; i<relation_num; i++)
        {
            for (int ii=0; ii<n; ii++)
                relation_vec[i][ii] = randn(0,1.0/n,-1,1);
        }
        for (int i=0; i<entity_num; i++)
        {
            for (int ii=0; ii<n; ii++)
                entity_vec[i][ii] = randn(0,1.0/n,-1,1);
        }

        bfgs();
    }

private:
    int n;
	int method;
    double res;//loss function value
    double count,count1;//loss function gradient
    double rate;//learning rate
    double belta;
	double margin;
    vector<int> fb_h,fb_l,fb_r;
    vector<vector<int> > feature;
    vector<vector<double> > A, A_tmp;
    vector<vector<double> > relation_vec,entity_vec,relation_tmp,entity_tmp;
    double norm(vector<double> &a)
    {
        double x = vec_len(a);
        if (x>1)
        for (int ii=0; ii<a.size(); ii++)
                a[ii]/=x;
        return 0;
    }
    double norm2one(vector<double> &a)
    {
        double x = vec_len(a);
        for (int ii=0; ii<a.size(); ii++)
                a[ii]/=x;
        return 0;
    }
    double norm(vector<double> &a, vector<double> &A)
    {
		norm2one(A);
    	double sum=0;
    	while (true)
    	{
			for (int i=0; i<n; i++)
				sum+=sqr(A[i]);
			sum = sqrt(sum);
			for (int i=0; i<n; i++)
				A[i]/=sum;
		    double x=0;
		    for (int ii=0; ii<n; ii++)
		    {
		        x+=A[ii]*a[ii];
		    }
		    if (x>0.1)
		    {
		        for (int ii=0; ii<n; ii++)
		        {
		        	a[ii]-=rate*A[ii];
					A[ii]-=rate*a[ii];
		        }
		    }
		    else
		    	break;
		}
		norm2one(A);
        return 0;
    }
    int rand_max(int x)
    {
        int res = (rand()*rand())%x;
        if (res<0)
            res+=x;
        return res;
    }

    void bfgs()
    {
        int step=0;
        count=0;
        count1=0;
        res=0;
        double times=0;
        int nbatches=50;
        int neval = 1000;
        int batchsize = fb_h.size()/nbatches;

        A_tmp = A;
        relation_tmp = relation_vec;
        entity_tmp = entity_vec;
        for (int eval=0; eval<neval; eval++)
        {

        	res=0;
         	for (int batch = 0; batch<nbatches; batch++)
         	{
         		for (int k=0; k<batchsize; k++)
         		{
					int i=rand_max(fb_h.size());
					int j=rand_max(entity_num);
					double pr = 1000*right_mean[fb_r[i]]/(right_mean[fb_r[i]]+left_mean[fb_r[i]]);
					if (method ==0)
                        pr = 500;
					if (rand()%1000<pr)
					{
						while (ok[make_pair(fb_h[i],fb_r[i])].count(j)>0)
							j=rand_max(entity_num);
                        train_kb(fb_h[i],fb_l[i],fb_r[i],fb_h[i],j,fb_r[i]);
					}
					else
					{
						while (ok[make_pair(j,fb_r[i])].count(fb_l[i])>0)
							j=rand_max(entity_num);
                        train_kb(fb_h[i],fb_l[i],fb_r[i],j,fb_l[i],fb_r[i]);
					}
                	norm(entity_tmp[fb_h[i]]);
                	norm(entity_tmp[fb_l[i]]);
                	norm(entity_tmp[j]);
					norm(entity_tmp[fb_h[i]],A_tmp[fb_r[i]]);
					norm(entity_tmp[fb_l[i]],A_tmp[fb_r[i]]);
					norm(entity_tmp[j],A_tmp[fb_r[i]]);
         		}
				
           		A = A_tmp;
           		relation_vec = relation_tmp;
           		entity_vec = entity_tmp;
           	}
           		cout<<eval<<' '<<res<<endl;
                FILE* f1 = fopen(("A.txt"+version).c_str(),"w");
                FILE* f2 = fopen(("relation2vec.txt"+version).c_str(),"w");
                FILE* f3 = fopen(("entity2vec.txt"+version).c_str(),"w");
                times+=1;
                for (int i=0; i<relation_num; i++)
                {
                    for (int ii=0; ii<n; ii++)
                        fprintf(f2,"%.6lf\t",relation_vec[i][ii]);
                    fprintf(f2,"\n");
                }
                for (int i=0; i<relation_num; i++)
                {
                	for (int jj=0; jj<n; jj++)
                		fprintf(f1,"%.6lf\t",A[i][jj]);
                	fprintf(f1,"\n");
                }
                for (int i=0; i<entity_num; i++)
                {
                    for (int ii=0; ii<n; ii++)
                        fprintf(f3,"%.6lf\t",entity_vec[i][ii]);
                    fprintf(f3,"\n");
                }
                fclose(f1);
                fclose(f2);
                fclose(f3);
        }
    }
    double res1;
    double calc_sum(int e1,int e2,int rel)
    {
        double tmp1=0,tmp2=0;
        for (int jj=0; jj<n; jj++)
        {
        	tmp1+=A[rel][jj]*entity_vec[e1][jj];
            tmp2+=A[rel][jj]*entity_vec[e2][jj];
        }

        double sum=0;
        for (int ii=0; ii<n; ii++)
            sum+=fabs(entity_vec[e2][ii]-tmp2*A[rel][ii]-(entity_vec[e1][ii]-tmp1*A[rel][ii])-relation_vec[rel][ii]);
        return sum;
    }
    void gradient(int e1,int e2,int rel,double belta)
    {
        double tmp1 = 0, tmp2 = 0;
        double sum_x=0;
        for (int jj=0; jj<n; jj++)
        {
            tmp1+=A[rel][jj]*entity_vec[e1][jj];
            tmp2+=A[rel][jj]*entity_vec[e2][jj];
        }
        for (int ii=0; ii<n; ii++)
        {

            double x = 2*(entity_vec[e2][ii]-tmp2*A[rel][ii]-(entity_vec[e1][ii]-tmp1*A[rel][ii])-relation_vec[rel][ii]);
			//for L1 distance function
            if (x>0)
            	x=1;
            else
            	x=-1;
            sum_x+=x*A[rel][ii];
            relation_tmp[rel][ii]-=belta*rate*x;
            entity_tmp[e1][ii]-=belta*rate*x;
            entity_tmp[e2][ii]+=belta*rate*x;
            A_tmp[rel][ii]+=belta*rate*x*tmp1;
            A_tmp[rel][ii]-=belta*rate*x*tmp2;
        }
        for (int ii=0; ii<n; ii++)
        {
        	A_tmp[rel][ii]+=belta*rate*sum_x*entity_vec[e1][ii];
        	A_tmp[rel][ii]-=belta*rate*sum_x*entity_vec[e2][ii];
        }

        norm(relation_tmp[rel]);
		norm(entity_tmp[e1]);
		norm(entity_tmp[e2]);

        norm2one(A_tmp[rel]);
        norm(relation_tmp[rel],A_tmp[rel]);
    }
    void train_kb(int e1_a,int e2_a,int rel_a,int e1_b,int e2_b,int rel_b)
    {
        double sum1 = calc_sum(e1_a,e2_a,rel_a);
        double sum2 = calc_sum(e1_b,e2_b,rel_b);
        if (sum1+margin>sum2)
        {
        	res+=margin+sum1-sum2;
        	gradient( e1_a, e2_a, rel_a, -1);
        	gradient(e1_b, e2_b, rel_b,1);
        }
    }
};

Train train;
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
	while (fscanf(f_kb,"%s",buf)==1)
    {
        string s1=buf;
        fscanf(f_kb,"%s",buf);
        string s2=buf;
        fscanf(f_kb,"%s",buf);
        string s3=buf;
        if (entity2id.count(s1)==0)
        {
            cout<<"miss entity:"<<s1<<endl;
        }
        if (entity2id.count(s2)==0)
        {
            cout<<"miss entity:"<<s2<<endl;
        }
        if (relation2id.count(s3)==0)
        {
            relation2id[s3] = relation_num;
            relation_num++;
        }
        int e1 = entity2id[s1];
        int e2 = entity2id[s2];
        int rel = relation2id[s3];
        if (left_candidate_ok[rel].count(e1)==0)
        {
            left_candidate_ok[rel][e1]=1;
            left_candidate[rel].push_back(e1);
        }
        if (right_candidate_ok[rel].count(e2)==0)
        {
            right_candidate_ok[rel][e2]=1;
            right_candidate[rel].push_back(e2);
        }
        entity2num[e1]++;
        entity2num[e2]++;
        left_entity[rel][e1].push_back(e2);
        right_entity[rel][e2].push_back(e1);
        train.add(e1,e2,rel);
    }
    for (int i=0; i<relation_num; i++)
    {
    	double sum1=0,sum2=0,sum3 = 0;
    	for (map<int,vector<int> >::iterator it = left_entity[i].begin(); it!=left_entity[i].end(); it++)
    	{
    		sum1++;
    		sum2+=it->second.size();
    		sum3+=sqr(it->second.size());
    	}
    	left_mean[i]=sum2/sum1;

    	left_var[i]=sum3/sum1-sqr(left_mean[i]);
    }
    for (int i=0; i<relation_num; i++)
    {
    	double sum1=0,sum2=0,sum3=0;
    	for (map<int,vector<int> >::iterator it = right_entity[i].begin(); it!=right_entity[i].end(); it++)
    	{
    		sum1++;
    		sum2+=it->second.size();
    		sum3+=sqr(it->second.size());
    	}
    	right_mean[i]=sum2/sum1;
    	right_var[i]=sum3/sum1-sqr(right_mean[i]);
    }

    for (int i=0; i<relation_num; i++)
    	cout<<i<<'\t'<<id2relation[i]<<' '<<left_mean[i]<<' '<<right_mean[i]<<endl;

    fclose(f_kb);


    cout<<"relation_num="<<relation_num<<endl;
    cout<<"entity_num="<<entity_num<<endl;
}

int ArgPos(char *str, int argc, char **argv) {
  int a;
  for (a = 1; a < argc; a++) if (!strcmp(str, argv[a])) {
    if (a == argc - 1) {
      printf("Argument missing for %s\n", str);
      exit(1);
    }
    return a;
  }
  return -1;
}

int main(int argc,char**argv)
{
    srand((unsigned) time(NULL));
    int method = 1;
    int n = 100;
    double rate = 0.001;
    double margin = 1;
    int i;
    if ((i = ArgPos((char *)"-size", argc, argv)) > 0) n = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-margin", argc, argv)) > 0) margin = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-method", argc, argv)) > 0) method = atoi(argv[i + 1]);
    cout<<"size = "<<n<<endl;
    cout<<"learing rate = "<<rate<<endl;
    cout<<"margin = "<<margin<<endl;
    if (method)
        version = "bern";
    else
        version = "unif";
    cout<<"method = "<<version<<endl;
    prepare();
    train.run(n,rate,margin,method);
}


