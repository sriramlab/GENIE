#ifndef IO_H

#define IO_H

#include "std.h"
#include "printable.h"

class io {
	public:
		static int debug ;
		static void print (const string& s, int level, string file) ; 
		static void println (const string& s, int level) ;
		static void print (const string& s, int level )  ;
		static void print (printable &p, int level) ;
		static void println(const printable &p, int level) ;
		static void print (printable &p, int level, string file) ;
};



inline
void io::print (const string &s, int level )  {
	if (level <= debug)
		cout << s << flush;
}

inline
void io::println (const string &s, int level) {
	print (s+"\n",level);
}

inline
void io::print (const string &s, int level, string file) { 
	if (level <= debug) { 
		ofstream ofs (file.c_str());
		ofs << s;
		ofs.close();
	}
}

inline
void io::print (printable &p, int level) { 
	print (p.to_string(), level);
}

inline
void io::println (const printable &p, int level) {
	io::print (p.to_string()+"\n", level);
}

inline
void io::print (printable &p, int level, string file) { 
	print (p.to_string(), level, file);
}



template<class Matrix>
void write_matrix(std::ofstream &ofs, const Matrix& matrix){
	int rows = static_cast<int>(matrix.rows());
    int cols = static_cast<int>(matrix.cols());
    ofs.write((char*) (&rows), sizeof(int));
    ofs.write((char*) (&cols), sizeof(int));
    ofs.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
}

template<class Matrix>
void read_matrix(std::ifstream &ifs, Matrix& matrix){
	int rows = 0; int cols = 0;
    ifs.read((char*) (&rows),sizeof(int));
    ifs.read((char*) (&cols),sizeof(int));
    matrix.resize(rows, cols);
    ifs.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
}


#endif
