#include "std.h"
#include "functions.h"

bool functions::fileExists(const std::string& fileName)
{
		std::fstream fin;
		fin.open(fileName.c_str(),std::ios::in);
		if( fin.is_open() )
		{
				fin.close();
				return true;
		}
		fin.close();
		return false;
}


void functions::tokenize(const string& str,
				vector<string>& tokens,
				const string& delimiters)
{
		// Skip delimiters at beginning.
		string::size_type lastPos = str.find_first_not_of(delimiters, 0);
		//         // Find first "non-delimiter".
		string::size_type pos     = str.find_first_of(delimiters, lastPos);

		while (string::npos != pos || string::npos != lastPos)
		{
				// Found a token, add it to the vector.
				tokens.push_back(str.substr(lastPos, pos - lastPos));
				// Skip delimiters.  Note the "not_of"
				lastPos = str.find_first_not_of(delimiters, pos);
				// Find next "non-delimiter"
				pos = str.find_first_of(delimiters, lastPos);
		}
}




