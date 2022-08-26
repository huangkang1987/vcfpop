/* String Functions */

#include "vcfpop.h"

/*Read bcf genotype string to an array */
TARGET void ReadBCFGenoString(char* gt, ushort* alleles, bool& phased, int& v, int asize, int vlen, uint64 l, char* name)
{
	phased = true;
	for (; v < vlen; ++v)
	{
		switch (asize)
		{
		case 1:
			alleles[v] = *(byte*)gt++;
			if (alleles[v] == 0x80 || alleles[v] == 0x81) //padding for lower ploidy level
				return;
			break;
		case 2:
			alleles[v] = *(ushort*)gt; gt += 2;
			if (alleles[v] == 0x8000 || alleles[v] == 0x8001) //padding for lower ploidy level
				return;
			break;
		}
		if (v && (alleles[v] & 0x1) == 0) 
			phased = false;

		alleles[v] = (alleles[v] >> 1) - 1;

		if (name && alleles[v] >= locus[l].k)
			Exit("\nError: allele index exceed number of alleles (%d), in individual %s, variant %s\n", locus[l].k, name, locus[l].GetName());

		if (alleles[v] == 0xFFFF)
		{
			SetFF(alleles, N_MAX_PLOIDY);
			return;
		}
	}
}

/*Read vcf genotype string to an array */
TARGET void ReadVCFGenoString(ushort* alleles, char* genostr, int ploidy, uint64 l, char* name)
{
	for (int i = 0; i < ploidy; ++i)
	{
		int tallele = -1;
		if (*genostr != '.')
		{
			tallele = ReadInteger(genostr);
			if (name && tallele >= locus[l].k)
				Exit("\nError: allele index exceed number of alleles (%d), in individual %s, variant %s\n", locus[l].k, name, locus[l].GetName());
			genostr++;
		}

		if (tallele == -1)
		{
			SetFF(alleles, ploidy);
			break;
		}
		alleles[i] = (ushort)tallele;
	}
}

/* Append val to str */
TARGET void AppendString(char*& str, const char* _val)
{
	char* val = (char*)_val;
	while (*val) *str++ = *val++;
}

/* Read bcf typed int from string */
TARGET uint ReadTypedInt(char*& str)
{
	byte type = *(byte*)str++;
	uint re = 0;
	switch (type & 0xF)
	{
	case 1: re = *(byte*)str; str += 1; break;
	case 2: re = *(ushort*)str; str += 2; break;
	case 3: re = *(uint*)str; str += 4; break;
	default: Exit("\nError: can not read typed integer.\n");
	}
	return re;
}

/* write bcf typed int into string */
TARGET void AppendTypedInt(char*& str, uint val)
{
	if (val >> 8 == 0)
	{
		*str++ = 0x11;
		*str++ = (byte)val;
	}
	else if (val >> 16 == 0)
	{
		*str++ = 0x12;
		*(ushort*)str = (ushort)val;
		str += 2;
	}
	else
	{
		*str++ = 0x13;
		*(uint*)str = val;
		str += 4;
	}
}

/* Fast print an integer to a string */
TARGET void AppendInt(char*& str, int val)
{
	if (val < 0)
	{
		*str++ = '-';
		val = -val;
	}
	int mask = 1000000000;
	while (mask > val) mask /= 10;
	if (!mask)
	{
		*str++ = '0';
		return;
	}
	while (mask)
	{
		*str++ = (char)('0' + (val / mask));
		val %= mask;
		mask /= 10;
	}
}

/* Print with indent */
TARGET void printi(const char* a)
{
	int64 len = strlen(a);
	char abuf[121];
	while (len > 0)
	{
		printf("    ");
		int64 plen = Min(len, 116ll);
		memmove(abuf, a, plen);
		abuf[plen] = '\0';
		printf((const char*)abuf);
		len -= 116;
		a += 116;
	}
}

/* Trim double quotes in a string and after -par= */
TARGET string TrimQuote(const string& a)
{
	if (a.size() == 0) return a;

	string re = a;
	if (re[0] == '\'' || re[0] == '\"')
		re = re.substr(1, re.size() - 1);

	if (re[re.size() - 1] == '\'' || re[re.size() - 1] == '\"')
		re = re.substr(0, re.size() - 1);

	int idx = (int)re.find_first_of('=');
	if (idx != -1 && re.size() > idx && (re[idx + 1] == '\'' || re[idx + 1] == '\"'))
		re.erase(idx + 1, 1);

	return re;
}

/* Trim double quotes and remove -par= in a string */
TARGET string TrimParQuote(const string& a)
{
	string re = a.substr(a.find_first_of('=') + 1);;

	if (re[0] == '\'' || re[0] == '\"')
		re = re.substr(1, re.size() - 1);

	if (re[re.size() - 1] == '\'' || re[re.size() - 1] == '\"')
		re = re.substr(0, re.size() - 1);

	return re;
}

/* Parse range '1-2' */
TARGET void ParseTwoNumber(char* names, uint& v1, uint& v2)
{
	v1 = v2 = 0xFFFFFFFF;
	if (names[0] == '#')
		sscanf(names, "#%d-#%d", &v1, &v2);
}

/* Parse range '1-2' */
TARGET void ParseTwoNumber(char* names, int& v1, int& v2)
{
	v1 = v2 = -1;
	if (names[0] == '#')
		sscanf(names, "#%d-#%d", &v1, &v2);
}

/* Compare two lines, and termined if A ends with linebreak comma | \0 or specified deliminator */
TARGET int LineCmpAterm(const char* a, const char* b, char termin)
{
	for (; ; a++, b++)
	{
		if (*a == termin || *a == '\r' || *a == '\n' || *a == '|' || *a == ',' || !*a)
			return 0;
		if (*a != *b)
			return 1;
	}
}

/* Compare two lines, and termined if A or B end with linebreak comma | \0 or specified deliminator */
TARGET int LineCmpABterm(const char* a, const char* b, char termin)
{
	for (; ; a++, b++)
	{
		if ((*a == termin || *a == '\r' || *a == '\n' || *a == '|' || *a == ',' || !*a) &&
			(*b == termin || *b == '\r' || *b == '\n' || *b == '|' || *b == ',' || !*b))
			return 0;
		if (*a != *b)
			return 1;
	}
}

/* Compare two lines, and termined if A ends with linebreak  |  \0 */
TARGET int LineCmp(const char* a, const char* b)
{
	for (; ; a++, b++)
	{
		if (*a == '\r' || *a == '\n' || *a == '|' || *a == ',' || !*a)
			return 0;
		if (*a != *b)
			return 1;
	}
}

/* Compare the lower case of two lines, and termined if A ends with linebreak or \0 */
TARGET int LwrLineCmp(const char* a, const char* b)
{
	for (; ; a++, b++)
	{
		if (*a == '\r' || *a == '\n' || !*a)
			return 0;
		if (*a + ((*a >= 65 && *a <= 90) ? 32 : 0) !=
			*b + ((*b >= 65 && *b <= 90) ? 32 : 0))
			return 1;
	}
}

/* Compare the lower case of two lines, and termined if A ends with linebreak or \0 */
TARGET int LwrLineCmp(const string& a, const string& b)
{
	for (int i = 0; ; ++i)
	{
		if (a[i] == '\r' || a[i] == '\n' || !a[i])
			return 0;
		if (a[i] + ((a[i] >= 65 && a[i] <= 90) ? 32 : 0) !=
			b[i] + ((b[i] >= 65 && b[i] <= 90) ? 32 : 0))
			return 1;
	}
}

/* Compare the lower case of two parameters */
TARGET int LwrParCmp(const char* a, const char* b)
{
	for (; ; a++, b++)
	{
		if ((*a == '\r' || *a == '\n' || !*a || *a == '|' || *a == ',') ||
			(*b == '\r' || *b == '\n' || !*b || *b == '|' || *b == ','))
			return !((*a == '\r' || *a == '\n' || !*a || *a == '|' || *a == ',') &&
				(*b == '\r' || *b == '\n' || !*b || *b == '|' || *b == ','));
		if (*a + ((*a >= 65 && *a <= 90) ? 32 : 0) != *b + ((*b >= 65 && *b <= 90) ? 32 : 0))
			return 1;
	}
}

/* Compare the lower case of two string */
TARGET int LwrStrCmp(const char* a, const char* b)
{
	for (; ; a++, b++)
	{
		if (!*a || !*b)
			return !(!*a && !*b);
		if (*a + ((*a >= 65 && *a <= 90) ? 32 : 0) !=
			*b + ((*b >= 65 && *b <= 90) ? 32 : 0))
			return 1;
	}
}

/* Compare the lower case of two string */
TARGET int LwrStrCmp(const string& a, const string& b)
{
	for (int i = 0; ; ++i)
	{
		if (!a[i] || !b[i])
			return !(!a[i] && !b[i]);
		if (a[i] + ((a[i] >= 65 && a[i] <= 90) ? 32 : 0) !=
			b[i] + ((b[i] >= 65 && b[i] <= 90) ? 32 : 0))
			return 1;
	}
}

/* Read all text of a file */
TARGET string ReadAllText(const string& file)
{
	FILE* f1 = fopen(file.c_str(), "rb");
	if (!f1) Exit("\nError: Cannot open file %s.\n", file.c_str());
	fseeko64(f1, 0ll, SEEK_END);
	uint64 flen = ftello64(f1);
	char* buf = new char[flen + 1];
	fseeko64(f1, 0ll, SEEK_SET);
	fread(buf, 1, flen, f1);
	fclose(f1);

	buf[flen] = '\0';
	string re = buf;
	delete[] buf;
	return re;
}

/* Replace substrings */
TARGET string ReplaceStr(const string& str, const string& a, const string& b, char escape)
{
	string re = "";
	int p1 = 0, alen = (int)a.size();
	bool isescape = true;
	while (str[p1])
	{
		if (str[p1] == escape) isescape = !isescape;
		if (isescape && str[p1] == a[0] && !memcmp(a.c_str(), str.c_str() + p1, alen))
		{
			p1 += alen;
			re += b;
		}
		else
			re += str[p1++];
	}
	return re;
}

/* Replace characters, not allocate memory */
TARGET void ReplaceChar(char* str, char a, char b)
{
	for (;;)
	{
		if (*str == a) *str = b;
		if (*++str == '\0') return;
	}
}

/* Does the string has a char */
TARGET bool ContainsChar(char* str, char ch, int len)
{
	for (int i = 0; i < len; ++i)
		if (str[i] == ch)
			return true;
	return false;
}

/* Does the string has a char */
TARGET bool ContainsChar(char* str, char ch)
{
	for (; *str; str++)
		if (*str == ch)
			return true;
	return false;
}

/* Count char in a string */
TARGET int CountChar(char* str, char ch)
{
	int count = 0;
	for (; *str; str++)
		count += (*str == ch);
	return count;
}

/* Count char val in a string */
TARGET int CountChar(const string& str, char ch)
{
	return (int)std::count(str.begin(), str.end(), ch);
}

/* Count any elements of val in a string */
TARGET int CountChars(char* str, const char* ch)
{
	uint64 count = 0;
	for (; *str; str++)
		for (uint64 i = 0; ch[i]; ++i)
			if (*str == ch[i] && ++count) break;
	return (int)count;
}

/* Count any elements of val in a string */
TARGET int CountChars(char* str, const char* ch, int64 len)
{
	uint64 count = 0;
	for (; len; --len, str++)
		for (uint64 i = 0; ch[i]; ++i)
			if (*str == ch[i] && ++count) break;
	return (int)count;
}

/* Skip rep lines */
TARGET char* LineNextIdx(char* str, const char* ch, int64 rep)
{
	char* val = (char*)ch;
	for (; ; str++)
	{
		if (*str == *val && !LwrLineCmp(val, str) && !--rep) return str;
		if (*str == '\0' || *str == '\n') return NULL;
	}
}

/* Skip rep vals */
TARGET char* StrNextIdx(char* str, const char* ch, int64 rep)
{
	char* val = (char*)ch;
	for (; ; str++)
	{
		if (*str == *val && !LwrLineCmp(val, str) && !--rep) return str;
		if (!*str) return NULL;
	}
}

/* Skip rep vals */
TARGET char* StrNextIdx(char* str, char ch, int64 rep)
{
	for (; ; str++)
	{
		if (*str == ch && !--rep) return str;
		if (!*str) return NULL;
	}
}

/* Next space char */
TARGET char* StrNextSpace(char* str)
{
	if (!*str) return NULL;

	//is space, return
	if (*str == ' ' || *str == '\t') return str;

	//next non-space char
	while (*str && (*str == ' ' || *str == '\t')) str++;
	if (!*str) return NULL;

	//next space
	while (*str && (*str != ' ' && *str != '\t')) str++;
	if (!*str) return NULL;

	return str;
}

/* Next non-space char */
TARGET char* StrNextChar(char* str)
{
	if (!*str) return NULL;

	//is non-space char, return
	if (*str != ' ' && *str != '\t') return str;

	//next space
	while (*str && (*str != ' ' && *str != '\t')) str++;
	if (!*str) return NULL;

	//next non-space char
	while (*str && (*str == ' ' || *str == '\t')) str++;
	if (!*str) return NULL;

	return str;
}

/* Skip rep \0s */
TARGET char* StrNextIdx0(char* str, int64 rep)
{
	for (; ; str++)
		if (*str == '\0' && !--rep) return str;
}

/* Replace delim with \0 and return pointers (new allocated) to each string */
TARGET char** SplitStr(char* str, char delim, int& count)
{
	count = (uint)CountChar(str, delim) + 1;
	char** re = new char*[count];
	SetZero(re, count);
	count = 0;
	bool inquote = false;
	while (*str)
	{
		while (!inquote && *str == delim) str++;
		if (*str == '\"')
		{
			inquote = true;
			memmove(str, str + 1, strlen(str));
		}

		if (*str) re[count++] = str;
		else break;
		while (*str && (inquote || *str != delim))
		{
			if (*str == '\"')
			{
				inquote = !inquote;
				memmove(str, str + 1, strlen(str));
				str--;
			}
			str++;
		}
		if (*str == delim) *str++ = '\0';
		else break;
	}
	return re;
}

/* Replace delim with \0 and return pointers (new allocated) to each string */
TARGET char** SplitStr(char* str, char delim, int64& count)
{
	count = CountChar(str, delim) + 1;
	char** re = new char*[count];
	SetZero(re, count);
	count = 0;
	bool inquote = false;
	while (*str)
	{
		while (!inquote && *str == delim) str++;
		if (*str == '\"')
		{
			inquote = true;
			memmove(str, str + 1, strlen(str));
		}

		if (*str) re[count++] = str;
		else break;
		while (*str && (inquote || *str != delim))
		{
			if (*str == '\"')
			{
				inquote = !inquote;
				memmove(str, str + 1, strlen(str));
				str--;
			}
			str++;
		}
		if (*str == delim) *str++ = '\0';
		else break;
	}
	return re;
}

/* Split a string */
TARGET vector<string> SplitStr(const string& str, const string& delim)
{
	size_t pos_start = 0, pos_end, delim_len = delim.length();
	string token;
	vector<string> res;

	while ((pos_end = str.find(delim, pos_start)) != string::npos) {
		token = str.substr(pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back(token);
	}

	res.push_back(str.substr(pos_start));
	return res;
}

/* Convert to lower case characters */
TARGET void StrLwr(char* str)
{
	do
	{
		if (*str >= 65 && *str <= 90)
			*str += 32;
	} while (*str++);
}

/* Read an allele from spagedi files*/
TARGET int ReadIntegerSpagedi(char*& str, int maxdigit)
{
	while ((*str < '0' || *str > '9') && *str != '-') str++;

	int sign = 0, re = 0;
	if (*str == '-')
	{
		str++;
		sign = -1;
	}
	else if (*str >= '0' && *str <= '9')
	{
		sign = 1;
	}
	else
	{
		str++;
		return -1;
	}

	int md = 0;
	while (*str >= '0' && *str <= '9')
	{
		re = re * 10 + (*str++ - '0');
		if (++md == maxdigit) break;
	}

	return re * sign;
}

/* Fast read an integer from a string */
TARGET int64 ReadLong(char*& str)
{
	while ((*str < '0' || *str > '9') && *str != '-') str++;

	int64 sign = 0, re = 0;
	if (*str == '-')
	{
		str++;
		sign = -1;
	}
	else if (*str >= '0' && *str <= '9')
	{
		sign = 1;
	}
	else
	{
		str++;
		return -1;
	}

	while (*str >= '0' && *str <= '9')
		re = re * 10ll + (*str++ - '0');
	return re * sign;
}

/* Parse an integer in binary format */
TARGET uint ReadBinInteger(char*& str, int len)
{
	uint re = 0;
	switch (len)
	{
	case 1: re = *(byte*)str; break;
	case 2: re = *(ushort*)str; break;
	case 3: re = *(uint*)str & 0xFFFFFF; break;
	case 4: re = *(uint*)str; break;
	}
	str += len;
	return re;
}

/* Parse a integer from string and move pointer */
TARGET int ReadInteger(char*& str)
{
	while ((*str < '0' || *str > '9') && *str != '-' && *str != '?') str++;
	
	if (*str == '?')
	{
		while (*str++ == '?');
		return -1;
	}

	int sign = 0, re = 0;
	if (*str == '-')
	{
		str++;
		sign = -1;
	}
	else if (*str >= '0' && *str <= '9')
	{
		sign = 1;
	}
	else
	{
		str++;
		return -1;
	}

	while (*str >= '0' && *str <= '9')
		re = re * 10 + (*str++ - '0');
	return re * sign;
}

/* Parse a integer from string and do not move pointer */
TARGET int ReadIntegerKeep(char* str)
{
	while ((*str < '0' || *str > '9') && *str != '-') str++;

	int sign = 0, re = 0;
	if (*str == '-')
	{
		str++;
		sign = -1;
	}
	else if (*str >= '0' && *str <= '9')
	{
		sign = 1;
	}
	else
	{
		str++;
		return -1;
	}

	while (*str >= '0' && *str <= '9')
		re = re * 10 + (*str++ - '0');
	return re * sign;
}

/* Parse a real number from string and move pointer */
TARGET double ReadDouble(char*& str)
{
	int sign = 0;
	double re = 0, mask = 1;
	if (*str == '-')
	{
		str++;
		sign = -1;
	}
	else if (*str >= '0' && *str <= '9')
		sign = 1;
	else
	{
		str++;
		return -1;
	}

	while (*str >= '0' && *str <= '9')
		re = re * 10 + (*str++ - '0');
	if (*str == '.')
	{
		str++;
		while (*str >= '0' && *str <= '9')
		{
			mask *= 0.1;
			re = re + (*str++ - '0') * mask;
		}
	}
	return re * sign;
}

/* Parse a real number from string and do not move pointer */
TARGET double ReadDoubleKeep(char* str)
{
	int sign = 0;
	double re = 0, mask = 1;
	if (*str == '-')
	{
		str++;
		sign = -1;
	}
	else if (*str >= '0' && *str <= '9')
	{
		sign = 1;
	}
	else
	{
		str++;
		return -1;
	}

	while (*str >= '0' && *str <= '9')
		re = re * 10 + (*str++ - '0');
	if (*str == '.')
	{
		str++;
		while (*str >= '0' && *str <= '9')
		{
			mask *= 0.1;
			re = re + (*str++ - '0') * mask;
		}
	}
	return re * sign;
}

/* Read real range parameter */
TARGET void GetRangeParDouble(string gpar, bool& parid, double& minv, double& maxv, double min, double max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e1];
	gpar[e1] = '\0';
	sscanf(gpar.c_str() + s1, "%lf", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e2];
	gpar[e2] = '\0';
	sscanf(gpar.c_str() + s2, "%lf", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %lf and %lf.\n", gpar.c_str(), min, max);
}

/* Read 32bit integer range parameter */
TARGET void GetRangeParInteger(string gpar, bool& parid, uint& minv, uint& maxv, uint min, uint max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e1];
	sscanf(gpar.c_str() + s1, "%u", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e2];
	sscanf(gpar.c_str() + s2, "%u", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %u and %u.\n", gpar.c_str(), min, max);
}

/* Read 32bit integer range parameter */
TARGET void GetRangeParInteger(string gpar, bool& parid, int& minv, int& maxv, int min, int max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e1];
	sscanf(gpar.c_str() + s1, "%d", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e2];
	sscanf(gpar.c_str() + s2, "%d", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %d and %d.\n", gpar.c_str(), min, max);
}

/* Read 64bit integer range parameter */
TARGET void GetRangeParLong(string gpar, bool& parid, uint64& minv, uint64& maxv, uint64 min, uint64 max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e1];
	sscanf(gpar.c_str() + s1, "%llu", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e2];
	sscanf(gpar.c_str() + s2, "%llu", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %llu and %llu.\n", gpar.c_str(), min, max);
}

/* Read 64bit integer range parameter */
TARGET void GetRangeParLong(string gpar, bool& parid, int64& minv, int64& maxv, int64 min, int64 max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e1];
	sscanf(gpar.c_str() + s1, "%lld", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	b = gpar[e2];
	sscanf(gpar.c_str() + s2, "%lld", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %lld and %lld.\n", gpar.c_str(), min, max);
}

/* Read real parameter */
TARGET void GetParDouble(string gpar, bool& parid, double& val, double min, double max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());

	sscanf(gpar.c_str() + s1, "%lf", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %lf and %lf.\n", gpar.c_str(), min, max);
}

/* Read 32bit integer parameter */
TARGET void GetParInteger(string gpar, bool& parid, uint& val, uint min, uint max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());

	sscanf(gpar.c_str() + s1, "%u", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %u and %u.\n", gpar.c_str(), min, max);
}

/* Read 32bit integer parameter */
TARGET void GetParInteger(string gpar, bool& parid, int& val, int min, int max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());

	sscanf(gpar.c_str() + s1, "%d", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %d and %d.\n", gpar.c_str(), min, max);
}

/* Read 64bit integer parameter */
TARGET void GetParLong(string gpar, bool& parid, uint64& val, uint64 min, uint64 max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	sscanf(gpar.c_str() + s1, "%llu", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %llu and %llu.\n", gpar.c_str(), min, max);
}

/* Read 64bit integer parameter */
TARGET void GetParLong(string gpar, bool& parid, int64& val, int64 min, int64 max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());
	sscanf(gpar.c_str() + s1, "%lld", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %lld and %lld.\n", gpar.c_str(), min, max);
}

/* Read boolean parameter */
TARGET void GetParBool(string gpar, bool& parid)
{
	if (parid)
		Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;
}

/* Read string parameter */
TARGET void GetParString(string gpar, const string& ref, bool& parid, int& val)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;

	int idx = (int)gpar.find_first_of('=');
	if (idx == -1)
		Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());

	gpar = gpar.substr(idx + 1);
	vector<string> refs = SplitStr(ref, "|");

	auto it = std::find(refs.begin(), refs.end(), gpar);
	if (it == refs.end())
		Exit("\nError: Unrecognized parameter: %s.\n ", gpar.c_str());
	val = std::distance(refs.begin(), it) + 1;
}

/* Read multiple section parameters */
TARGET void GetParStringMultiSel(string gpar, const string& ref, bool& parid, byte* val)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar.c_str());
	parid = true;

	int idx = (int)gpar.find_first_of('=');
	if (idx == -1)
		Exit("\nError: cannot parse parameter %s, check format.\n", gpar.c_str());

	vector<string> pars = SplitStr(gpar.substr(idx + 1), ",");
	vector<string> refs = SplitStr(ref, "|");
	memset(val, 0, N_MAX_OPTION);

	for (int i = 0; i < pars.size(); ++i)
	{
		auto it = std::find(refs.begin(), refs.end(), pars[i]);
		if (it == refs.end())
			Exit("\nError: Unrecognized option: %s in parameter: \n%s\n.", pars[i].c_str(), gpar.c_str());
		val[std::distance(refs.begin(), it) + 1] = 1;
	}
}

/* Print a real number to a file */
TARGET void WriteReal(FILE* fout, double val)
{
	if (IsError(val)) fprintf(fout, "nan");
	else fprintf(fout, g_decimal_str, val);
}

/* Print a real number to a string */
TARGET void WriteReal(char*& fout, double val)
{
	if (IsError(val)) sprintf(fout, "nan");
	else sprintf(fout, g_decimal_str, val);
	fout += strlen(fout);
}

/* Append a real number to a string */
TARGET void AppendReal(char* sout, double val)
{
	if (IsError(val)) sprintf(sout, "nan");
	else sprintf(sout, g_decimal_str, val);
}
