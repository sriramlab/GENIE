#ifndef VECTORFN_H

#define VECTORFN_H

namespace vectorfn {

template <class T>
inline void init(vector<T> &t, T x)  {
	for (int i = 0; i < t.size(); i++)
		t[i] = x;
}


template <class T> 
inline void add (vector<T> &t, T x) { 
	for (int i = 0; i < t.size(); i++)
		t[i] += x;

}

template <class T> 
inline void sub (vector<T> &t, T x) { 
	for (int i = 0; i < t.size(); i++)
		t[i] -= x;
}

template <class T> 
inline void mul (vector<T> &t, T x) { 
	for (int i = 0; i < t.size(); i++)
		t[i] *= x;
}

template <class T> 
inline void div (vector<T> &t, T x) { 
	if (x==0) {
		cerr << "Error: Divide by zero\n";
		exit(1);
	}
	for (int i = 0; i < t.size(); i++)
		t[i] /= x;

}

template <class T> 
inline T iprod (vector<T> &t, vector<T> &s) { 
	if (t.size() != s.size()) {
		cerr << "Incompatible sizes\n";
		exit(1);
	}
	T sum = 0 ;
	for (int i = 0; i < t.size(); i++)
		sum += t[i]*s[i];
	return sum;
}

template <class T> 
inline void idiv (vector<T> &t, vector<T> &x) { 
	for (int i = 0; i < t.size(); i++) {
		if (x[i] == 0) {
			cerr << "Error: Divide by zero\n";
			exit(1);
		}
		t[i] /= x[i];
	}
}

template <class T> 
inline T sum (vector<T> &t) { 
	T sum = 0 ;
	for (int i = 0; i < t.size(); i++)
		sum += t[i];
	return sum;

}

template <class T> 
inline T mean (vector<T> &t) { 
	T mean = 0 ;
	for (int i = 0; i < t.size(); i++)
		mean += t[i];
    if (t.size () >  0 )
        mean/= t.size();
	return mean;

}

template <class T> 
inline pair<T,int> min (vector<T> &t) { 
	T m  (0) ;
	int ind = -1;
	bool flag = false;
	for (int i = 0; i < t.size(); i++){
		if (!std::isinf(t[i]) && !std::isnan(t[i])){
			if (!flag) {
				m = t[i];
				ind = i;
				flag = true;
			} else if (t[i] < m) {
				m = t[i];
				ind = i;
			}
		}
	}
	return pair<T, int> (m,ind);

}

template <class T> 
inline pair<T,int> max (vector<T> &t) { 
	T m  (0) ;
	int ind = -1;
	bool flag = false;
	for (int i = 0; i < t.size(); i++){
		if (!std::isinf(t[i]) && !std::isnan(t[i])){
			if (!flag) {
				m = t[i];
				ind = i;
				flag = true;
			} else if (t[i] > m) {
				m = t[i];
				ind = i;
			}
		}
	}
	return pair<T, int> (m,ind);

}


template <class T> 
inline T lsumexp (vector<T> &t) {
	pair<T,int> p = max(t);
	T sum = 0;
	if (p.second>=0){
		for (int i = 0; i < t.size(); i++) {
			if (!std::isinf(t[i]) && !std::isnan(t[i]))
				sum += exp(t[i]-p.first);
		}
		sum = log(sum) + p.first;

	} else {
		if (t.size() >0 )
			return log(sum);
		else 
			cerr << "Empty vector"<<endl;
	}
	return sum;
}


template <class T>
inline void printvectornl(vector<T> &t, string delim = " "){
    printvector (t, delim, true);
}

template <class T>
inline void printvector(vector<T> &t, string delim = " ", bool newline = false){
		for (int i = 0; i < t.size(); i++)
				cout << t[i] << delim;
        if (newline)
            cout << endl;
}


/*
template <class T>
inline void quantile (vector<T> &t, vector<double> &q, vector<double> &out) { 
    vector<T> tmp (t.size());
    std::copy ( t.begin(), t.end(), tmp.begin());
    std::sort (tmp.begin(),tmp.end());
    out.resize (q.size());
    
    for (int j  = 0 ;  j < q.size(); j++) {
        int q1 = ceil(t.size()*q[j]);
        int q2 = floor (t.size()*q[j]+1);
        out[j] = 0.5*(tmp[q1-1]+tmp[q2-1]);
    }
}

}*/

inline void quantile (vector<double> &t, vector<double> &q, vector<double> &out) { 
    vector<double> tmp (t.size());
    std::copy ( t.begin(), t.end(), tmp.begin());
    std::sort (tmp.begin(),tmp.end());
    out.resize (q.size());
    
    for (int j  = 0 ;  j < q.size(); j++) {
        int q1 = ceil(t.size()*q[j]);
        int q2 = floor (t.size()*q[j]+1);
        if (q1 >= 1 && q2 >= 1 )
            out[j] = 0.5*(tmp[q1-1]+tmp[q2-1]);
        else 
            out[j] = 0;
    }
}

}


#endif
