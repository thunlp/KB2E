
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

bool L1_flag=1;

//normal distribution
double rand(double min, double max)
{
    return min+(max-min)*rand()/(RAND_MAX+1.0);
}
double normal(double x, double miu,double sigma)
{
    return 1.0/sqrt(2*pi)/sigma*exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}
double sigmod(double x)
{
	return 1.0/(1+exp(-2*x));
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
	res = sqrt(res);
	return res;
}
void vec_output(vector<double> a)
{
	for (int i=0; i<a.size(); i++)
		cout<<a[i]<<"\t";
	cout<<endl;
}

string version;
char buf[100000],buf1[100000];
int relation_num,entity_num;
map<string,int> relation2id,entity2id;
map<int,string> id2entity,id2relation;

map<int,int> entity2num;


map<int,map<int,vector<int> > > left_entity,right_entity;
map<int,double> left_mean,right_mean,left_var,right_var;

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
        m = n_in;
        rate = rate_in;
        margin = margin_in;
        method = method_in;
		A.resize(relation_num);
		for (int i=0; i<relation_num; i++)
		{
		    A[i].resize(n);
		    for (int jj=0; jj<n; jj++)
		    {
		        A[i][jj].resize(m);
		        for (int ii=0; ii<m; ii++)
		        {
					if (ii==jj)
						A[i][jj][ii]=1;
					else
						A[i][jj][ii]=0;
		        }
		    }
		}
        relation_vec.resize(relation_num);
		for (int i=0; i<relation_vec.size(); i++)
			relation_vec[i].resize(m);
        entity_vec.resize(entity_num);
		for (int i=0; i<entity_vec.size(); i++)
			entity_vec[i].resize(n);
        for (int i=0; i<relation_num; i++)
        {
            for (int ii=0; ii<m; ii++)
                relation_vec[i][ii] = randn(0,1.0/m,-1,1);
        }
        FILE* f1 = fopen("../TransE/entity2vec.unif","r");
        for (int i=0; i<entity_num; i++)
        {
            for (int ii=0; ii<n; ii++)
            	fscanf(f1,"%lf",&entity_vec[i][ii]);
            norm(entity_vec[i]);
        }
        fclose(f1);

		FILE* f2 = fopen("../TransE/relation2vec.unif","r");
        for (int i=0; i<relation_num; i++)
        {
            for (int ii=0; ii<n; ii++)
            	fscanf(f2,"%lf",&relation_vec[i][ii]);
        }
        fclose(f2);
        bfgs();
    }

private:
    int n,m,method;
    double res;//loss function value
    double rate,margin;//learning rate
    vector<int> fb_h,fb_l,fb_r;
    vector<vector<double> > relation_vec,entity_vec;
    vector<vector<double> > relation_tmp,entity_tmp;
    vector<vector<vector<double> > > A,A_tmp;
    void norm(vector<double> &a)
    {
        double x = vec_len(a);
        if (x>1)
        for (int ii=0; ii<a.size(); ii++)
                a[ii]/=x;
    }
    void norm(vector<double> &a, vector<vector<double> > &A)
    {
    	while (true)
    	{
		    double x=0;
		    for (int ii=0; ii<m; ii++)
		    {
		        double tmp = 0;
		        for (int jj=0; jj<n; jj++)
		            tmp+=A[jj][ii]*a[jj];
		        x+=sqr(tmp);
		    }
		    if (x>1)
		    {
		        double lambda=1;
		        for (int ii=0; ii<m; ii++)
		        {
		            double tmp = 0;
		            for (int jj=0; jj<n; jj++)
		                tmp+=A[jj][ii]*a[jj];
		            tmp*=2;
		            for (int jj=0; jj<n; jj++)
		            {
		                A[jj][ii]-=rate*lambda*tmp*a[jj];
		                a[jj]-=rate*lambda*tmp*A[jj][ii];
		            }
		        }
		    }
		    else
		    	break;
		}
    }
    int rand_max(int x)
    {
        int res = (rand()*rand())%x;
        while (res<0)
            res+=x;
        return res;
    }

    void bfgs()
    {
        res=0;
        int nbatches=100;
        int nepoch = 1000;
        int batchsize = fb_h.size()/nbatches;
        relation_tmp=relation_vec;
        entity_tmp = entity_vec;
        A_tmp = A;
            for (int epoch=0; epoch<nepoch; epoch++)
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
                		norm(relation_tmp[fb_r[i]]);
                		norm(entity_tmp[fb_h[i]]);
						norm(entity_tmp[fb_l[i]]);
						norm(entity_tmp[j]);
						norm(entity_tmp[fb_h[i]],A_tmp[fb_r[i]]);
						norm(entity_tmp[fb_l[i]],A_tmp[fb_r[i]]);
						norm(entity_tmp[j],A_tmp[fb_r[i]]);
						norm(entity_tmp[k]);
						norm(entity_tmp[k],A_tmp[fb_r[i]]);
             		}
		            relation_vec = relation_tmp;
		            entity_vec = entity_tmp;
		            A = A_tmp;
             	}
                cout<<"epoch:"<<epoch<<' '<<res<<endl;
                FILE* f1 = fopen(("relation2vec."+version).c_str(),"w");
                FILE* f2 = fopen(("entity2vec."+version).c_str(),"w");
                FILE* f3 = fopen(("A."+version).c_str(),"w");
                for (int i=0; i<relation_num; i++)
                {
                    for (int ii=0; ii<m; ii++)
                        fprintf(f1,"%.6lf\t",relation_vec[i][ii]);
                    fprintf(f1,"\n");
                }
                for (int i=0; i<entity_num; i++)
                {
                    for (int ii=0; ii<n; ii++)
                        fprintf(f2,"%.6lf\t",entity_vec[i][ii]);
                    fprintf(f2,"\n");
                }
                for (int i=0; i<relation_num; i++)
                for (int jj=0; jj<n; jj++)
                {
                    for (int ii=0; ii<m; ii++)
                    {
                        fprintf(f3,"%.6lf\t",A[i][jj][ii]);
                    }
                    fprintf(f3,"\n");
                }
                fclose(f1);
                fclose(f2);
                fclose(f3);
            }
    }
    double res1;
    double calc_sum(int e1,int e2,int rel,int same)
    {
    	vector<double> e1_vec;
        e1_vec.resize(m);
        vector<double> e2_vec;
        e2_vec.resize(m);
        for (int ii=0; ii<m; ii++)
        {
            for (int jj=0; jj<n; jj++)
            {
                e1_vec[ii]+=A[rel][jj][ii]*entity_vec[e1][jj];
                e2_vec[ii]+=A[rel][jj][ii]*entity_vec[e2][jj];
            }
        }
        double sum=0;
        if (L1_flag)
        	for (int ii=0; ii<m; ii++)
            	sum+=fabs(e2_vec[ii]-e1_vec[ii]-same*relation_vec[rel][ii]);
        else
        	for (int ii=0; ii<m; ii++)
            	sum+=sqr(e2_vec[ii]-e1_vec[ii]-same*relation_vec[rel][ii]);
        return sum;
    }
    void gradient_one(int e1, int e2, int rel, int belta,int same)
    {
    	for (int ii=0; ii<m; ii++)
        {
            double tmp1 = 0, tmp2 = 0;
            for (int jj=0; jj<n; jj++)
            {
                tmp1+=A[rel][jj][ii]*entity_vec[e1][jj];
                tmp2+=A[rel][jj][ii]*entity_vec[e2][jj];
            }
            double x = 2*(tmp2-tmp1-relation_vec[rel][ii]);
            if (L1_flag)
            	if (x>0)
            		x=1;
            	else
            		x=-1;
            for (int jj=0; jj<n; jj++)
            {
                A_tmp[rel][jj][ii]-=belta*rate*x*(entity_vec[e1][jj]-entity_vec[e2][jj]);
                entity_tmp[e1][jj]-=belta*rate*x*A[rel][jj][ii];
                entity_tmp[e2][jj]+=belta*rate*x*A[rel][jj][ii];
            }
            relation_tmp[rel][ii]-=same*belta*rate*x;
        }
    }
    void gradient(int e1_a,int e2_a,int rel_a,int e1_b,int e2_b,int rel_b)
    {
    	gradient_one(e1_a,e2_a,rel_a,-1,1);
    	gradient_one(e1_b,e2_b,rel_b,1,1);
    }
    void train_kb(int e1_a,int e2_a,int rel_a,int e1_b,int e2_b,int rel_b)
    {
        double sum1 = calc_sum(e1_a,e2_a,rel_a,1);
        double sum2 = calc_sum(e1_b,e2_b,rel_b,1);
        if (sum1+margin>sum2)
        {
        	res+=margin+sum1-sum2;
        	gradient( e1_a, e2_a, rel_a, e1_b, e2_b, rel_b);
        }
    }
    void gradient_same(int e1, int e2, int rel)
    {
    	double ll = entity2num[e2]*1.0/(entity2num[e1]+entity2num[e2]);
    	double rr = entity2num[e1]*1.0/(entity2num[e1]+entity2num[e2]);
    	for (int ii=0; ii<m; ii++)
        {
            double tmp1 = 0, tmp2 = 0;
            for (int jj=0; jj<n; jj++)
            {
                tmp1+=A[rel][jj][ii]*entity_vec[e1][jj];
                tmp2+=A[rel][jj][ii]*entity_vec[e2][jj];
            }
            double x = 2*(tmp2-tmp1-relation_vec[rel][ii]);
            if (L1_flag)
            	if (x>0)
            		x=1;
            	else
            		x=-1;
            for (int jj=0; jj<n; jj++)
            {
                A_tmp[rel][jj][ii]-=-1*rate*x*(ll*entity_vec[e1][jj]-rr*entity_vec[e2][jj]);
                entity_tmp[e1][jj]-=-1*rate*ll*x*A[rel][jj][ii];
                entity_tmp[e2][jj]+=-1*rate*rr*x*A[rel][jj][ii];
            }
        }
    }
    void same_kb(int e1,int e2,int rel)
    {
    	double sum = calc_sum(e1,e2,rel,0);
    	if (sum>1)
    	{
    		res+=sum-1;
		 	gradient_same(e1,e2,rel);
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
        entity2num[entity2id[s1]]++;
        entity2num[entity2id[s2]]++;
        left_entity[relation2id[s3]][entity2id[s1]].push_back(entity2id[s2]);
        right_entity[relation2id[s3]][entity2id[s2]].push_back(entity2id[s1]);
        train.add(entity2id[s1],entity2id[s2],relation2id[s3]);
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
    int n = 50;
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


