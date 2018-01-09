#include <bits/stdc++.h>
using namespace std;

#define PI acos(-1.0)
#define INF 1e30
#define eps 1e-12

int channel_span, total_states;
vector<double> prior, cluster_mean;
vector<vector<double>> transition_prob;

struct channel{
	vector<double> h;
	double noise_mean, noise_variance;

	default_random_engine generator;
	normal_distribution<double> distribution;

	void set(){
		distribution = normal_distribution<double>(noise_mean, noise_variance);
	}

	double transmit(string seq){
		double ret = 0;
		for(int i=seq.size()-1; i>=0; i--){
			if(seq[i] == '0') ret -= h[h.size()-i-1];
			else ret += h[h.size()-i-1];
		}
  		return ret + distribution(generator); 
	}

	void display(){
		puts("Weights: ");
		for(auto v: h) cout << v << ' ';
		cout << endl;

		cout << "noise_mean: " << noise_mean << "\nnoise_variance: " << noise_variance << "\n\n";
	}
};

void readConfig(string fileName, channel &ch)
{
	int N;
	double w;
	FILE *fp = fopen(fileName.c_str(), "r");

	fscanf(fp, "%lf %lf", &ch.noise_mean, &ch.noise_variance);
	ch.set();

	while(fscanf(fp, "%lf", &w) == 1){
		ch.h.push_back(w);
	}

	fclose(fp);
}

int getStateNumber(string state)
{
	int ret = 0;
	for(int i=0; i<state.size(); i++){
		ret = 2*ret + state[i]-'0';
	}
	return ret;
}

string getState(int state_number)
{
	string ret = "";
	while(state_number > 0){
		int dig = state_number%2;
		state_number /= 2;
		ret.push_back(dig+'0');
	}

	while(ret.size() < channel_span) ret.push_back('0');
	reverse(ret.begin(), ret.end());
	return ret;
}

string readFile(string fileName)
{
	ifstream fin;
	fin.open(fileName);

	string ret;
	getline(fin, ret);
	fin.close();

	return ret;
}

void normalize(vector<double> &v)
{
	double sum = 0;
	for(auto val: v) sum += val;
	for(auto &val: v) val /= sum;
}

void setParameters(string bseq, channel ch)
{
	///calculate prior probability
	prior.resize(total_states, 0);
	for(int i=0; i<bseq.size()-channel_span+1; i++){
		string state = "";
		for(int j=0; j<channel_span; j++) state.push_back(bseq[i+j]);
		int state_number = getStateNumber(state);
		prior[state_number]++;
	}

	normalize(prior);	

	//calculate transition probability
	transition_prob.resize(total_states, vector<double>(total_states, 0));

	string cur_state = "", next_state = "";
	int cur_state_number, next_state_number;

	for(int i=0; i<channel_span; i++) cur_state.push_back(bseq[i]);
	cur_state_number = getStateNumber(cur_state);

	for(int i=1; i<bseq.size()-channel_span+1; i++){
		for(int j=0; j<channel_span; j++) next_state.push_back(bseq[i+j]);

		next_state_number = getStateNumber(next_state);
		transition_prob[cur_state_number][next_state_number]+=1;
		
		cur_state = next_state;
		cur_state_number = next_state_number;
		next_state = "";
	}

	for(auto &v: transition_prob) normalize(v);

	///calculate cluster means
	vector<int> cluster_count(total_states, 0);
	cluster_mean.resize(total_states, 0);
	for(int i=0; i<bseq.size()-channel_span+1; i++){
		string state = "";
		for(int j=0; j<channel_span; j++) state.push_back(bseq[i+j]);

		int state_number = getStateNumber(state);
		double channel_output = ch.transmit(state);

		cluster_count[state_number]++;
		cluster_mean[state_number] += channel_output;
	}

	for(int i=0; i<total_states; i++) cluster_mean[i] /= cluster_count[i];

	/*
	for(auto v: cluster_mean) cout << v << ' ';
	cout << endl;
	*/
}

string viterbi(string bseq, channel ch)
{
	vector<vector<double>> D(total_states, vector<double>(bseq.size()-channel_span+1, 0));
	vector<vector<int>> savePath(total_states, vector<int>(bseq.size()-channel_span+1, -1));

	int N = bseq.size()-channel_span+1;

	string state = "";
	for(int i=0; i<channel_span; i++) state.push_back(bseq[i]);
	double x = ch.transmit(state);

	for(int i=0; i<total_states; i++){
		double val = x - cluster_mean[i];
		D[i][0] = log(prior[i]) - 0.5 * log(2*PI*ch.noise_variance) - val * val / (2.0 * ch.noise_variance);
	}

	for(int i=1; i<N; i++){
		string state = "";
		for(int j=0; j<channel_span; j++) state.push_back(bseq[i+j]);
		double x = ch.transmit(state);

		for(int j=0; j<total_states; j++){
			double val = x - cluster_mean[j];
			double mxval = -INF;

			for(int k=0; k<total_states; k++){
				if(abs(transition_prob[k][j]) < eps) continue; //transition not possible

				double cand = D[k][i-1] + log(transition_prob[k][j]) - 0.5 * log(2*PI*ch.noise_variance) - 
								val * val /(2.0 * ch.noise_variance);

				if(cand > mxval){
					savePath[j][i] = k;
					mxval = cand;
				}
			}

			D[j][i] = mxval;
		}
	}
    /*
    for(auto v: D){
        for(auto u: v){
            cout << u << ' ';
        }
        cout << endl;
    }
    */

	int last_class = 0;
	for(int i=1; i<total_states; i++){
		if(D[i][N-1] > D[last_class][N-1]){
			last_class = i;
		}
	}

	vector<int> path;
	path.push_back(last_class);

	for(int i=N-1; i>0; i--){
		path.push_back(savePath[last_class][i]);
		last_class = savePath[last_class][i];
	}

	
	//cout << "path size: " << path.size() << endl;
	/*
	for(auto v: path) cout << v << ' ';
	cout << endl;
	*/

	string ans = getState(path.back());
	for(int i=path.size()-2; i>=0; i--) ans.push_back(getState(path[i]).back());
	return ans;
}

int main()
{
	channel ch;
	readConfig("config.txt", ch);
	//ch.display();
	
	channel_span = ch.h.size();
	total_states = (int)pow(2, channel_span);
	//cout << total_states << endl;
	
    string bseq = readFile("train.txt");
	setParameters(bseq, ch);

	string target = readFile("test.txt");
	string original = viterbi(target, ch);

	//cout << target.size() << endl;
	//cout << original.size() << endl;

	if(target == original){
		cout << "Complete Match.\n";
	}
	else{
		cout << "Error exists.\n";
	}

	FILE *fp = fopen("out.txt", "w");
	fprintf(fp, "%s", original.c_str());
	fclose(fp);

	//cout << original << endl;
	
	return 0;
}
