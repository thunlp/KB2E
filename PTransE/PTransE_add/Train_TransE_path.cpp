#include<iostream>
#include<cstring>
#include<cstdio>
#include<map>
#include<vector>
#include<string>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include<sstream>
using namespace std;


#define pi 3.1415926535897932384626433832795


map<vector<int>,string> path2s;


map<pair<string,int>,double>  path_confidence;

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
/*	if (L1_flag)
		for (int i=0; i<a.size(); i++)
			res+=fabs(a[i]);
	else*/
	{
		for (int i=0; i<a.size(); i++)
			res+=a[i]*a[i];
		res = sqrt(res);
	}
	return res;
}

string version;
char buf[100000],buf1[100000],buf2[100000];
int relation_num,entity_num;
map<string,int> relation2id,entity2id;
map<int,string> id2entity,id2relation;




vector<vector<pair<int,int> > > path;

class Train{

public:
	map<pair<int,int>, map<int,int> > ok;
    void add(int x,int y,int z, vector<pair<vector<int>,double> > path_list)
    {
        fb_h.push_back(x);
        fb_r.push_back(z);
        fb_l.push_back(y);
		fb_path.push_back(path_list);
        ok[make_pair(x,z)][y]=1;
    }
    void run()
    {
        n = 100;
        rate = 0.001;
		cout<<"n="<<n<<' '<<"rate="<<rate<<endl;
		relation_vec.resize(relation_num);
		for (int i=0; i<relation_vec.size(); i++)
			relation_vec[i].resize(n);
        entity_vec.resize(entity_num);
		for (int i=0; i<entity_vec.size(); i++)
			entity_vec[i].resize(n);
        relation_tmp.resize(relation_num);
		for (int i=0; i<relation_tmp.size(); i++)
			relation_tmp[i].resize(n);
        entity_tmp.resize(entity_num);
		for (int i=0; i<entity_tmp.size(); i++)
			entity_tmp[i].resize(n);
        for (int i=0; i<relation_num; i++)
        {
            for (int ii=0; ii<n; ii++)
                relation_vec[i][ii] = randn(0,1.0/n,-6/sqrt(n),6/sqrt(n));
        }
        for (int i=0; i<entity_num; i++)
        {
            for (int ii=0; ii<n; ii++)
                entity_vec[i][ii] = randn(0,1.0/n,-6/sqrt(n),6/sqrt(n));
            norm(entity_vec[i]);
        }

        bfgs();
    }

private:
    int n;
    double res;//loss function value
    double count,count1;//loss function gradient
    double rate;//learning rate
    double belta;
    vector<int> fb_h,fb_l,fb_r;
	vector<vector<pair<vector<int>,double> > >fb_path;
    vector<vector<int> > feature;
    vector<vector<double> > relation_vec,entity_vec;
    vector<vector<double> > relation_tmp,entity_tmp;
	
    double norm(vector<double> &a)
    {
        double x = vec_len(a);
        if (x>1)
        for (int ii=0; ii<a.size(); ii++)
                a[ii]/=x;
        return 0;
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
        double margin = 1,margin_rel = 1;
        cout<<"margin="<<' '<<margin<<"margin_rel="<<margin_rel<<endl;
        res=0;
        int nbatches=100;
        int neval = 1000;
        int batchsize = fb_h.size()/nbatches;
 		relation_tmp=relation_vec;
		entity_tmp = entity_vec;
		map<string, int> fb_count;
        for (int eval=0; eval<neval; eval++)
        {
        	res=0;
         	for (int batch = 0; batch<nbatches; batch++)
         	{
				int e1 = rand_max(entity_num);
         		for (int k=0; k<batchsize; k++)
         		{

					int j=rand_max(entity_num);
					int i=rand_max(fb_h.size());
					int e1 = fb_h[i], rel = fb_r[i], e2  = fb_l[i];
					
					int rand_tmp = rand()%100;
					if (rand_tmp<25)
					{
						while (ok[make_pair(e1,rel)].count(j)>0)
							j=rand_max(entity_num);
                        train_kb(e1,e2,rel,e1,j,rel,margin);
					}
					else
					if (rand_tmp<50)
					{
						while (ok[make_pair(j,rel)].count(e2)>0)
							j=rand_max(entity_num);
                        train_kb(e1,e2,rel,j,e2,rel,margin);

					}
					else
					{
						int rel_neg = rand_max(relation_num);
						while (ok[make_pair(e1,rel_neg)].count(e2)>0)
							rel_neg = rand_max(relation_num);
                        train_kb(e1,e2,rel,e1,e2,rel_neg,margin);
					}
					if (fb_path[i].size()>0)
					{
						int rel_neg = rand_max(relation_num);
						while (ok[make_pair(e1,rel_neg)].count(e2)>0)
							rel_neg = rand_max(relation_num);
						for (int path_id = 0; path_id<fb_path[i].size(); path_id++)
						{
							vector<int> rel_path = fb_path[i][path_id].first;
							string  s = "";
							if (path2s.count(rel_path)==0)
							{
							    ostringstream oss;//创建一个流
								for (int ii=0; ii<rel_path.size(); ii++)
								{
									oss<<rel_path[ii]<<" ";
								}
							    s=oss.str();//
								path2s[rel_path] = s;
							}
							s = path2s[rel_path];

							double pr = fb_path[i][path_id].second;
							double pr_path = 0;
							if (path_confidence.count(make_pair(s,rel))>0)
								pr_path = path_confidence[make_pair(s,rel)];
							pr_path = 0.99*pr_path + 0.01;
							train_path(rel,rel_neg,rel_path,2*margin,pr*pr_path);
						}
					}
					norm(relation_tmp[rel]);
            		norm(entity_tmp[e1]);
            		norm(entity_tmp[e2]);
            		norm(entity_tmp[j]);
					e1 = e2;
         		}
	            relation_vec = relation_tmp;
	            entity_vec = entity_tmp;
         	}
            cout<<"eval:"<<eval<<' '<<res<<endl;
            FILE* f2 = fopen(("relation2vec.txt"+version).c_str(),"w");
            FILE* f3 = fopen(("entity2vec.txt"+version).c_str(),"w");
            for (int i=0; i<relation_num; i++)
            {
                for (int ii=0; ii<n; ii++)
                    fprintf(f2,"%.6lf\t",relation_vec[i][ii]);
                fprintf(f2,"\n");
            }
            for (int i=0; i<entity_num; i++)
            {
                for (int ii=0; ii<n; ii++)
                    fprintf(f3,"%.6lf\t",entity_vec[i][ii]);
                fprintf(f3,"\n");
            }
            fclose(f2);
            fclose(f3);
        }
    }
    double res1;
    double calc_kb(int e1,int e2,int rel)
    {
        double sum=0;
        for (int ii=0; ii<n; ii++)
		{
			double tmp = entity_vec[e2][ii]-entity_vec[e1][ii]-relation_vec[rel][ii];
	        if (L1_flag)
				sum+=fabs(tmp);
			else
				sum+=sqr(tmp);
		}
        return sum;
    }
    void gradient_kb(int e1,int e2,int rel, double belta)
    {
        for (int ii=0; ii<n; ii++)
        {

            double x = 2*(entity_vec[e2][ii]-entity_vec[e1][ii]-relation_vec[rel][ii]);
            if (L1_flag)
            	if (x>0)
            		x=1;
            	else
            		x=-1;
            relation_tmp[rel][ii]-=belta*rate*x;
            entity_tmp[e1][ii]-=belta*rate*x;
            entity_tmp[e2][ii]+=belta*rate*x;
        }
    }
    double calc_path(int r1,vector<int> rel_path)
    {
        double sum=0;
        for (int ii=0; ii<n; ii++)
		{
			double tmp = relation_vec[r1][ii];
			for (int j=0; j<rel_path.size(); j++)
				tmp-=relation_vec[rel_path[j]][ii];
	        if (L1_flag)
				sum+=fabs(tmp);
			else
				sum+=sqr(tmp);
		}
        return sum;
    }
    void gradient_path(int r1,vector<int> rel_path, double belta)
    {
        for (int ii=0; ii<n; ii++)
        {

			double x = relation_vec[r1][ii];
			for (int j=0; j<rel_path.size(); j++)
				x-=relation_vec[rel_path[j]][ii];
            if (L1_flag)
            	if (x>0)
            		x=1;
            	else
            		x=-1;
            relation_tmp[r1][ii]+=belta*rate*x;
			for (int j=0; j<rel_path.size(); j++)
            	relation_tmp[rel_path[j]][ii]-=belta*rate*x;
        }
    }
    void train_kb(int e1_a,int e2_a,int rel_a,int e1_b,int e2_b,int rel_b,double margin)
    {
        double sum1 = calc_kb(e1_a,e2_a,rel_a);
        double sum2 = calc_kb(e1_b,e2_b,rel_b);
        if (sum1+margin>sum2)
        {
        	res+=margin+sum1-sum2;
        	gradient_kb(e1_a, e2_a, rel_a, -1);
			gradient_kb(e1_b, e2_b, rel_b, 1);
        }
    }
    void train_path(int rel, int rel_neg, vector<int> rel_path, double margin,double x)
    {
        double sum1 = calc_path(rel,rel_path);
        double sum2 = calc_path(rel_neg,rel_path);
		double lambda = 1;
        if (sum1+margin>sum2)
        {
        	res+=x*lambda*(margin+sum1-sum2);
        	gradient_path(rel,rel_path, -x*lambda);
			gradient_path(rel_neg,rel_path, x*lambda);
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
		id2relation[x+1345] = "-"+st;
		relation_num++;
	}
    FILE* f_kb = fopen("../data/train_pra.txt","r");
	while (fscanf(f_kb,"%s",buf)==1)
    {
        string s1=buf;
        fscanf(f_kb,"%s",buf);
        string s2=buf;
        if (entity2id.count(s1)==0)
        {
            cout<<"miss entity:"<<s1<<endl;
        }
        if (entity2id.count(s2)==0)
        {
            cout<<"miss entity:"<<s2<<endl;
        }
        int e1 = entity2id[s1];
        int e2 = entity2id[s2];
        int rel;
		fscanf(f_kb,"%d",&rel);
		fscanf(f_kb,"%d",&x);
		vector<pair<vector<int>,double> > b;
		b.clear();
		for (int i = 0; i<x; i++)
		{
			int y,z;
			vector<int> rel_path;
			rel_path.clear();
			fscanf(f_kb,"%d",&y);
			for (int j=0; j<y; j++)
			{
				fscanf(f_kb,"%d",&z);
				rel_path.push_back(z);
			}
			double pr;
			fscanf(f_kb,"%lf",&pr);
			b.push_back(make_pair(rel_path,pr));
		}
		//cout<<e1<<' '<<e2<<' '<<rel<<' '<<b.size()<<endl;
        train.add(e1,e2,rel,b);
    }
	relation_num*=2;
   
    cout<<"relation_num="<<relation_num<<endl;
    cout<<"entity_num="<<entity_num<<endl;
	
	FILE* f_confidence = fopen("../data/confidence.txt","r");
	while (fscanf(f_confidence,"%d",&x)==1)
	{
		string s = "";
		for (int i=0; i<x; i++)
		{
			fscanf(f_confidence,"%s",buf);
			s = s + string(buf)+" ";
		}
		fscanf(f_confidence,"%d",&x);
		for (int i=0; i<x; i++)
		{
			int y;
			double pr;
			fscanf(f_confidence,"%d%lf",&y,&pr);
		//	cout<<s<<' '<<y<<' '<<pr<<endl;
			path_confidence[make_pair(s,y)] = pr;
		}
	}
	fclose(f_confidence);
    fclose(f_kb);
}

int main(int argc,char**argv)
{
    srand((unsigned) time(NULL));
    if (argc!=2)
        return 0;
    else
    {
        version = argv[1];
        prepare();
        train.run();
    }
}


