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
int relation_num,entity_num,fb_relation_num;

map<int,int> cluster2fb;

map<string,int> relation2id,entity2id;
map<int,string> id2entity,id2relation;

map<int,int> entity2num;


map<int,map<int,int> > left_entity,right_entity;
map<int,double> left_mean,right_mean;

map<int,map<int,int> > fb_left_entity,fb_right_entity;
map<int,double> fb_left_mean,fb_right_mean;

class Train{

public:
	map<pair<int,int>, map<int,int> > ok,ok_fb;
    void add(int x,int y,int z)
    {
        fb_h.push_back(x);
        fb_r.push_back(z);
        fb_l.push_back(y);
        ok[make_pair(x,z)][y]=1;
        ok_fb[make_pair(x,cluster2fb[z])][y]=1;
    }
    double cluster_rate;
    void run(int n_in,double rate_in,double margin_in,int method_in)
    {
        n = n_in;
        m = n_in;
        rate = rate_in;
        margin = margin_in;
        method = method_in;
		A.resize(fb_relation_num);
		for (int i=0; i<fb_relation_num; i++)
		{
		    A[i].resize(n);
		    for (int jj=0; jj<n; jj++)
		        A[i][jj].resize(m);
		}

        relation_vec.resize(relation_num);
		for (int i=0; i<relation_vec.size(); i++)
        {
			relation_vec[i].resize(m);
			for (int ii=0; ii<m; ii++)
                relation_vec[i][ii] = randn(0,1.0/m,-1,1);
        }
        fb_relation_vec.resize(fb_relation_num);
		for (int i=0; i<fb_relation_vec.size(); i++)
        {
			fb_relation_vec[i].resize(m);
			for (int ii=0; ii<m; ii++)
                fb_relation_vec[i][ii] = randn(0,1.0/m,-1,1);
        }
        entity_vec.resize(entity_num);
		for (int i=0; i<entity_vec.size(); i++)
			entity_vec[i].resize(n);
        FILE* f1 = fopen("../TransW/entity2vec.unif","r");
        for (int i=0; i<entity_num; i++)
        {
            for (int ii=0; ii<n; ii++)
            	fscanf(f1,"%lf",&entity_vec[i][ii]);
            norm(entity_vec[i]);
        }
        fclose(f1);

		FILE* f2 = fopen("../TransW/relation2vec.unif","r");
        for (int i=0; i<fb_relation_num; i++)
        {
            for (int ii=0; ii<n; ii++)
            	fscanf(f2,"%lf",&fb_relation_vec[i][ii]);
        }
        fclose(f2);
        for (int i=0; i<relation_num; i++)
        {
            for (int ii=0; ii<n; ii++)
                relation_vec[i][ii] = fb_relation_vec[cluster2fb[i]][ii];
        }
        FILE* f3 = fopen("../TransW/A.unif","r");
        for (int i=0; i<fb_relation_num; i++)
            for (int jj=0; jj<m; jj++)
                for (int ii=0; ii<n; ii++)
                    fscanf(f3,"%lf",&A[i][jj][ii]);
        fclose(f3);
        bfgs();
    }
private:
    int n,m,method;
    double res;//loss function value
    double rate,margin;//learning rate
    vector<int> fb_h,fb_l,fb_r;
    vector<vector<double> > relation_vec,entity_vec, fb_relation_vec,relation_tmp,entity_tmp, fb_relation_tmp;
    vector<vector<vector<double> > > A,A_tmp, fb_A, fb_A_tmp;
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
		     //   res+=lambda*(x-1);
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
        int neval = 1000;
        int batchsize = fb_h.size()/nbatches;
        relation_tmp=relation_vec;
        entity_tmp = entity_vec;
        A_tmp = A;
        fb_relation_tmp = fb_relation_vec;
        fb_A_tmp = fb_A;

        for (int eval=0; eval<neval; eval++)
        {

            res=0;
            for (int batch = 0; batch<nbatches; batch++)
            {
                //cout<<batch<<endl;
                for (int kk=0; kk<batchsize; kk++)
                {
                    int i=rand_max(fb_h.size());
                    int j=rand_max(entity_num);
                    int e1 = fb_h[i];
                    int e2 = fb_l[i];
                    int rel = fb_r[i];
                    int fb_rel = cluster2fb[rel];
                    double pr = 1000*right_mean[fb_r[i]]/(right_mean[fb_r[i]]+left_mean[fb_r[i]]);
                    if (method ==0)
                        pr = 500;
                    if (rand()%1000<pr)
                    {
                        while (ok[make_pair(e1,rel)].count(j)>0||ok_fb[make_pair(e1,fb_rel)].count(j)>0)
                            j=rand_max(entity_num);
                        train_kb(e1,fb_l[i],rel,e1,j,rel);
                    }
                    else
                    {
                        while (ok[make_pair(j,rel)].count(fb_l[i])>0||ok_fb[make_pair(j,fb_rel)].count(fb_l[i])>0)
                            j=rand_max(entity_num);
                            train_kb(e1,fb_l[i],rel,j,fb_l[i],rel);
                    }
                    same_rel(rel);
                    norm(relation_tmp[rel]);
                    norm(fb_relation_tmp[fb_rel]);
                    norm(entity_tmp[e1]);
                    norm(entity_tmp[fb_l[i]]);
                    norm(entity_tmp[j]);
                    norm(entity_tmp[e1],A_tmp[fb_rel]);
                    norm(entity_tmp[fb_l[i]],A_tmp[fb_rel]);
                    norm(entity_tmp[j],A_tmp[fb_rel]);

                }
                relation_vec = relation_tmp;
                entity_vec = entity_tmp;
                A = A_tmp;
                fb_relation_tmp = fb_relation_vec;
                fb_A_tmp = fb_A;
            }
            cout<<"eval:"<<eval<<' '<<res<<endl;
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
            for (int i=0; i<fb_relation_num; i++)
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


            FILE* f4 = fopen(("fb_relation2vec."+version).c_str(),"w");
            for (int i=0; i<fb_relation_num; i++)
            {
                for (int ii=0; ii<m; ii++)
                    fprintf(f4,"%.6lf\t",fb_relation_vec[i][ii]);
                fprintf(f4,"\n");
            }
            fclose(f4);
        }
    }
    double res1;
    double calc_sum(int e1,int e2,int rel)
    {
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
            	sum+=fabs(e2_vec[ii]-e1_vec[ii]-relation_vec[rel][ii]*cluster_rate-(1-cluster_rate)*fb_relation_vec[fb_rel][ii]);
        else
        	for (int ii=0; ii<m; ii++)
            	sum+=sqr(e2_vec[ii]-e1_vec[ii]-relation_vec[rel][ii]*cluster_rate-(1-cluster_rate)*fb_relation_vec[fb_rel][ii]);
        return sum;
    }
    void gradient_one(int e1, int e2, int rel, int belta)
    {
        int fb_rel = cluster2fb[rel];
    	for (int ii=0; ii<m; ii++)
        {
            double tmp1 = 0, tmp2 = 0;
            for (int jj=0; jj<n; jj++)
            {
                tmp1+=A[fb_rel][jj][ii]*entity_vec[e1][jj];
                tmp2+=A[fb_rel][jj][ii]*entity_vec[e2][jj];
            }
            double x = 2*(tmp2-tmp1-relation_vec[rel][ii]*cluster_rate-(1-cluster_rate)*fb_relation_vec[fb_rel][ii]);
            if (L1_flag)
            	if (x>0)
            		x=1;
            	else
            		x=-1;
            for (int jj=0; jj<n; jj++)
            {
                A_tmp[fb_rel][jj][ii]-=belta*rate*x*(entity_vec[e1][jj]-entity_vec[e2][jj]);
                entity_tmp[e1][jj]-=belta*rate*x*A[fb_rel][jj][ii];
                entity_tmp[e2][jj]+=belta*rate*x*A[fb_rel][jj][ii];
            }
            relation_tmp[rel][ii]-=belta*rate*x*cluster_rate;
            fb_relation_tmp[fb_rel][ii]-=belta*rate*x*(1-cluster_rate);
        }
    }
    void gradient(int e1_a,int e2_a,int rel_a,int e1_b,int e2_b,int rel_b)
    {
    	gradient_one(e1_a,e2_a,rel_a,-1);
    	gradient_one(e1_b,e2_b,rel_b,1);
    }
    void train_kb(int e1_a,int e2_a,int rel_a,int e1_b,int e2_b,int rel_b)
    {
        double sum1 = calc_sum(e1_a,e2_a,rel_a);
        double sum2 = calc_sum(e1_b,e2_b,rel_b);
        if (sum1+margin>sum2)
        {
        	res+=margin+sum1-sum2;
        	gradient( e1_a, e2_a, rel_a, e1_b, e2_b, rel_b);
        }
    }
    void same_rel(double rel)
    {
        double sum=0;
        for (int ii=0; ii<n; ii++)
            if (L1_flag)
                sum+=fabs(fb_relation_vec[cluster2fb[rel]][ii]-relation_vec[rel][ii]);
            else
                sum+=sqr(fb_relation_vec[cluster2fb[rel]][ii]-relation_vec[rel][ii]);
        if (sum>0)
        {
            double lambda = 0.001;
            res+=lambda*(sum);
            for (int ii=0; ii<n; ii++)
            {
                double x = 2*(fb_relation_vec[cluster2fb[rel]][ii]-relation_vec[rel][ii]);
                if (L1_flag)
                    if (x>0)
                        x=1;
                    else
                        x=-1;
                fb_relation_tmp[cluster2fb[rel]][ii]-=lambda*rate*x;
                relation_tmp[rel][ii]+=lambda*rate*x;
            }
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
		fb_relation_num++;
	}
    FILE* f_kb = fopen("../data/train.txt","r");
    FILE* f_train_cluster = fopen("../cluster/output/train_cluster.txt","r");
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
            cout<<"miss rel:"<<s3<<endl;
        }
        int e1 = entity2id[s1];
        int e2 = entity2id[s2];
        int rel_fb = relation2id[s3];
        int rel;
        fscanf(f_train_cluster,"%d",&rel);
        left_entity[rel][e1]+=1;
        right_entity[rel][e2]+=1;
        if (cluster2fb.count(rel)==0)
        {
            relation_num+=1;
            cluster2fb[rel]=rel_fb;
        }
        else
        if (cluster2fb[rel]!=rel_fb)
        {
            cout<<"not map:"<<e1<<" "<<e2<<" "<<rel<<' '<<rel_fb<<endl;
            exit(-1);
        }

        fb_left_entity[cluster2fb[rel]][e1]+=1;
        fb_right_entity[cluster2fb[rel]][e2]+=1;
        train.add(e1,e2,rel);
    }
    fclose(f_train_cluster);
    fclose(f_kb);
    for (int i=0; i<relation_num; i++)
    {
    	double sum1=0,sum2=0;
    	for (map<int,int>::iterator it = left_entity[i].begin(); it!=left_entity[i].end(); it++)
    	{
    		sum1++;
    		sum2+=it->second;
    	}
    	left_mean[i]=sum2/sum1;
    }
    for (int i=0; i<relation_num; i++)
    {
    	double sum1=0,sum2=0;
    	for (map<int,int>::iterator it = right_entity[i].begin(); it!=right_entity[i].end(); it++)
    	{
    		sum1++;
    		sum2+=it->second;
    	}
    	right_mean[i]=sum2/sum1;
    }


     for (int i=0; i<fb_relation_num; i++)
    {
    	double sum1=0,sum2=0;
    	for (map<int,int>::iterator it = fb_left_entity[i].begin(); it!=fb_left_entity[i].end(); it++)
    	{
    		sum1++;
    		sum2+=it->second;
    	}
    	fb_left_mean[i]=sum2/sum1;
    }
    for (int i=0; i<fb_relation_num; i++)
    {
    	double sum1=0,sum2=0;
    	for (map<int,int>::iterator it = fb_right_entity[i].begin(); it!=fb_right_entity[i].end(); it++)
    	{
    		sum1++;
    		sum2+=it->second;
    	}
    	fb_right_mean[i]=sum2/sum1;
    }
    for (int i=0; i<fb_relation_num; i++)
    	cout<<fb_left_mean[i]<<' '<<fb_right_mean[i]<<endl;
    cout<<"fb_relation_num="<<fb_relation_num<<endl;
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


