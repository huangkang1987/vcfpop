/* String Functions */

#pragma once
#include "vcfpop.h"

/*Read bcf genotype string to an array */
TARGET void ReadBCFGenoString(char* gt, ushort* alleles, bool& phased, int& v, int asize, int vlen, uint64 l, char* name);

/*Read vcf genotype string to an array */
TARGET void ReadVCFGenoString(ushort* alleles, char* genostr, int ploidy, uint64 l, char* name);

/* Append val to str */
TARGET void AppendString(char*& str, const char* _val);

/* Read bcf typed int from string */
TARGET uint ReadTypedInt(char*& str);

/* write bcf typed int into string */
TARGET void AppendTypedInt(char*& str, uint val);

/* Fast print an integer to a string */
TARGET void AppendInt(char*& str, int val);

/* Print with indent */
TARGET void printi(const char* a);

/* Trim double quotes in a string and after -par= */
TARGET string TrimQuote(const string& a);

/* Trim double quotes and remove -par= in a string */
TARGET string TrimParQuote(const string& a);

/* Parse range '1-2' */
TARGET void ParseTwoNumber(char* names, uint& v1, uint& v2);

/* Parse range '1-2' */
TARGET void ParseTwoNumber(char* names, int& v1, int& v2);

/* Compare two lines, and termined if A ends with linebreak comma | \0 or specified deliminator */
TARGET int LineCmpAterm(const char* a, const char* b, char termin);

/* Compare two lines, and termined if A or B end with linebreak comma | \0 or specified deliminator */
TARGET int LineCmpABterm(const char* a, const char* b, char termin);

/* Compare two lines, and termined if A ends with linebreak  |  \0 */
TARGET int LineCmp(const char* a, const char* b);

/* Compare the lower case of two lines, and termined if A ends with linebreak or \0 */
TARGET int LwrLineCmp(const char* a, const char* b);

/* Compare the lower case of two lines, and termined if A ends with linebreak or \0 */
TARGET int LwrLineCmp(const string& a, const string& b);

/* Compare the lower case of two parameters */
TARGET int LwrParCmp(const char* a, const char* b);

/* Compare the lower case of two string */
TARGET int LwrStrCmp(const char* a, const char* b);

/* Compare the lower case of two string */
TARGET int LwrStrCmp(const string& a, const string& b);

/* Read all text of a file */
TARGET string ReadAllText(const string& file);

/* Replace substrings */
TARGET string ReplaceStr(const string& str, const string& a, const string& b, char escape = '\0');

/* Replace characters, not allocate memory */
TARGET void ReplaceChar(char* str, char a, char b);

/* Does the string has a char */
TARGET bool ContainsChar(char* str, char ch, int len);

/* Does the string has a char */
TARGET bool ContainsChar(char* str, char ch);

/* Count char in a string */
TARGET int CountChar(char* str, char ch);

/* Count char val in a string */
TARGET int CountChar(const string& str, char ch);

/* Count any elements of val in a string */
TARGET int CountChars(char* str, const char* ch);

/* Count any elements of val in a string */
TARGET int CountChars(char* str, const char* ch, int64 len);

/* Skip rep lines */
TARGET char* LineNextIdx(char* str, const char* ch, int64 rep);

/* Skip rep vals */
TARGET char* StrNextIdx(char* str, const char* ch, int64 rep);

/* Skip rep vals */
TARGET char* StrNextIdx(char* str, char ch, int64 rep);

/* Next space char */
TARGET char* StrNextSpace(char* str);

/* Next non-space char */
TARGET char* StrNextChar(char* str);

/* Skip rep \0s */
TARGET char* StrNextIdx0(char* str, int64 rep);

/* Replace delim with \0 and return pointers (new allocated) to each string */
TARGET char** SplitStr(char* str, char delim, int& count);

/* Replace delim with \0 and return pointers (new allocated) to each string */
TARGET char** SplitStr(char* str, char delim, int64& count);

/* Split a string */
TARGET vector<string> SplitStr(const string& str, const string& delim);

/* Convert to lower case characters */
TARGET void StrLwr(char* str);

/* Convert to lower case characters */
TARGET string StrLwr(string& str);

/* Read an allele from spagedi files*/
TARGET int ReadIntegerSpagedi(char*& str, int maxdigit);

/* Fast read an integer from a string */
TARGET int64 ReadLong(char*& str);

/* Parse an integer in binary format */
TARGET uint ReadBinInteger(char*& str, int len);

/* Parse a integer from string and move pointer */
TARGET int ReadInteger(char*& str);

/* Parse a integer from string and do not move pointer */
TARGET int ReadIntegerKeep(char* str);

/* Parse a real number from string and move pointer */
TARGET double ReadDouble(char*& str);

/* Parse a real number from string and do not move pointer */
TARGET double ReadDoubleKeep(char* str);

/* Read real range parameter */
TARGET void GetRangeParDouble(string gpar, bool& parid, double& minv, double& maxv, double min, double max);

/* Read 32bit integer range parameter */
TARGET void GetRangeParInteger(string gpar, bool& parid, uint& minv, uint& maxv, uint min, uint max);

/* Read 32bit integer range parameter */
TARGET void GetRangeParInteger(string gpar, bool& parid, int& minv, int& maxv, int min, int max);

/* Read 64bit integer range parameter */
TARGET void GetRangeParLong(string gpar, bool& parid, uint64& minv, uint64& maxv, uint64 min, uint64 max);

/* Read 64bit integer range parameter */
TARGET void GetRangeParLong(string gpar, bool& parid, int64& minv, int64& maxv, int64 min, int64 max);

/* Read real parameter */
TARGET void GetParDouble(string gpar, bool& parid, double& val, double min, double max);

/* Read 32bit integer parameter */
TARGET void GetParInteger(string gpar, bool& parid, uint& val, uint min, uint max);

/* Read 32bit integer parameter */
TARGET void GetParInteger(string gpar, bool& parid, int& val, int min, int max);

/* Read 64bit integer parameter */
TARGET void GetParLong(string gpar, bool& parid, uint64& val, uint64 min, uint64 max);

/* Read 64bit integer parameter */
TARGET void GetParLong(string gpar, bool& parid, int64& val, int64 min, int64 max);

/* Read 32bit integer array parameter */
TARGET void GetParIntegerArray(string gpar, bool& parid, uint* val, uint min, uint max, int npar);

/* Read 32bit integer array parameter */
TARGET void GetParIntegerArray(string gpar, bool& parid, int* val, int min, int max, int npar);

/* Read boolean parameter */
TARGET void GetParBool(string gpar, bool& parid);

/* Read string parameter */
TARGET void GetParString(string gpar, const string& ref, bool& parid, int& val);

/* Read string array parameter */
TARGET void GetParStringArray(string gpar, const string& ref, bool& parid, int* val, int npar);

/* Read multiple section parameters */
TARGET void GetParStringMultiSel(string gpar, const string& ref, bool& parid, byte* val);

/* Write a real number to a file */
TARGET void WriteReal(FILE* fout, double val);

/* Write a real number to a string */
TARGET void WriteReal(char*& fout, double val);

/* Append a real number to a string */
TARGET void AppendReal(char* sout, double val);

/* Write a real number to a file */
TARGET void WriteReal(FILE* fout, float val);

/* Write a real number to a string */
TARGET void WriteReal(char*& fout, float val);

/* Append a real number to a string */
TARGET void AppendReal(char* sout, float val);

/* Write a real number to a file in scientific notation */
TARGET void WriteScientific(FILE* fout, double val);

/* Write a real number to a string in scientific notation */
TARGET void WriteScientific(char*& fout, double val);

/* Append a real number to a string in scientific notation */
TARGET void AppendScientific(char* sout, double val);

/* Write a real number to a file in scientific notation */
TARGET void WriteScientific(FILE* fout, float val);

/* Write a real number to a string in scientific notation */
TARGET void WriteScientific(char*& fout, float val);

/* Append a real number to a string in scientific notation */
TARGET void AppendScientific(char* sout, float val);