/* File IO */

#include "vcfpop.h"

/* Map a file into memory */
TARGET byte* FileMapping::MapingReadOnlyFile(char* filename)
{
#ifdef _WIN64
	hFile = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);

	if (hFile == INVALID_HANDLE_VALUE)
		Exit("fail to map file %s into memory.", filename);

	DWORD size_low, size_high;
	size_low = GetFileSize(hFile, &size_high);

	hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, size_high, size_low, NULL);

	if (hMapFile == INVALID_HANDLE_VALUE)
		Exit("fail to map file %s into memory.", filename);

	addr = (byte*)MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0);

	return addr;
#else
	size = GetFileLen(filename);

	hFile = open(filename, O_RDONLY);

	if (hFile == -1)
		Exit("fail to map file %s into memory.", filename);

	addr = (byte*)mmap(NULL, size, PROT_READ, MAP_SHARED, hFile, 0);

	if (addr == MAP_FAILED)
		Exit("fail to map file %s into memory.", filename);

	return addr;
#endif
}

/* UnMap a file into memory */
TARGET void FileMapping::UnMapingReadOnlyFile()
{
#ifdef _WIN64
	UnmapViewOfFile(addr);
	CloseHandle(hMapFile);
	CloseHandle(hFile);
#else
	munmap(addr, size);
	close(hFile);
#endif
}

/* Get file size */
TARGET int64 GetFileLen(char* file)
{
	FILE* f1 = fopen(file, "rb");
	fseeko64(f1, 0, SEEK_END);
	int64 re = ftello64(f1);
	fclose(f1);
	return re;
}

/* Read a double from file */
TARGET double FGetDouble(FILE* file)
{
	double buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 8);
	else
		fread(&buf, 1, 8, file);
	return buf;
}

/* Read a float from file */
TARGET float FGetFloat(FILE* file)
{
	float buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 4);
	else
		fread(&buf, 1, 4, file);
	return buf;
}

/* Read an int64 from file */
TARGET int64 FGetInt64(FILE* file)
{
	int64 buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 8);
	else
		fread(&buf, 1, 8, file);
	return buf;
}

/* Read an uint64 from file */
TARGET uint64 FGetUint64(FILE* file)
{
	uint64 buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 8);
	else
		fread(&buf, 1, 8, file);
	return buf;
}

/* Read an int from file */
TARGET int FGetInt(FILE* file)
{
	int buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 4);
	else
		fread(&buf, 1, 4, file);
	return buf;
}

/* Read an uint from file */
TARGET uint FGetUint(FILE* file)
{
	uint buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 4);
	else
		fread(&buf, 1, 4, file);
	return buf;
}

/* Read a short from file */
TARGET short FGetShort(FILE* file)
{
	short buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 2);
	else
		fread(&buf, 1, 2, file);
	return buf;
}

/* Read an ushort from file */
TARGET ushort FGetUshort(FILE* file)
{
	ushort buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 2);
	else
		fread(&buf, 1, 2, file);
	return buf;
}

/* Read a byte from file */
TARGET byte FGetByte(FILE* file)
{
	byte buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 1);
	else
		fread(&buf, 1, 1, file);
	return buf;
}

/* Read a char from file */
TARGET char FGetChar(FILE* file)
{
	char buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 1);
	else
		fread(&buf, 1, 1, file);
	return buf;
}

/* Read a typed int from bcf file */
TARGET uint FGetTypedInt(FILE* file)
{
	byte type = FGetByte(file);
	switch (type & 0xF)
	{
	case 1: return FGetByte(file);
	case 2: return FGetUshort(file);
	case 3: return FGetUint(file);
	}
	return 0;
}

/* Is end of file */
TARGET bool FEof(FILE* file)
{
	if (g_format_val < 0)
		return gzeof((gzFile)file);
	else
		return feof(file);
}

/* Get compressed data offset */
TARGET int64 FOffset(FILE* file)
{
	if (g_format_val < 0)
		return gzoffset64((gzFile)file);
	else
		return ftello64(file);
}

/* Get decompressed data offset */
TARGET int64 FTell(FILE* file)
{
	if (g_format_val < 0)
		return gztell64((gzFile)file);
	else
		return ftello64(file);
}

/* Read decompressed data */
TARGET int64 FRead(void* buf, uint64 element, uint64 count, FILE* file)
{
	if (g_format_val < 0)
		return gzread((gzFile)file, buf, (uint)(element * count));
	else
		return fread(buf, element, count, file);
}

/* Change file offset for decompressed data */
TARGET int64 FSeek(FILE* file, int64 offset, int origin)
{
	if (g_format_val < 0)
		return gzseek64((gzFile)file, offset, origin);
	else
		return fseeko64(file, offset, origin);
}

/* Open an decompressed or compressed file */
TARGET FILE* FOpen(const char* file, const char* spec)
{
	if (g_format_val < 0)
	{
		gzFile f = gzopen(file, spec);
		gzbuffer(f, LINE_BUFFER);
		return (FILE*)f;
	}
	else
		return fopen(file, spec);
}

/* Close an decompressed or compressed file */
TARGET int FClose(FILE* file)
{
	if (g_format_val < 0)
		return gzclose((gzFile)file);
	else
		return fclose(file);
}

/* Format file name and open file */
TARGET FILE* FOpen(char* buf, const char* type, const char* fmt ...)
{
	va_list ap;
	int n = 0;
	va_start(ap, fmt);
	n = vsprintf(buf, fmt, ap);
	va_end(ap);
	FILE* f1 = fopen(buf, type);
	if (!f1) Exit("\nError: Cannot open file %s.\n", buf);
	return f1;
}

/* Open result file and write parameters */
TARGET void OpenResFile(const char* _spec, const char* title)
{
	char* spec = (char*)_spec;
	time_t time1;
	time(&time1);
	FRES_TIME = localtime(&time1);
	FRES = FOpen(FRES_NAME, "wb", "%s.%s.txt", g_output_val.c_str(), spec + 1);
	setvbuf(FRES, FRES_BUF, _IOFBF, LINE_BUFFER);
	const char* prefix = title[0] == '#' ? "#" : "";

	fprintf(FRES, "%s%s calculated by vcfpop v%s%s", prefix, title, VERSION, g_linebreak_val);
	fprintf(FRES, "%s  Input: %s%s", prefix, g_input_val.c_str(), g_linebreak_val);
	fprintf(FRES, "%s  Output: %s%s", prefix, g_output_val.c_str(), g_linebreak_val);
	fprintf(FRES, "%s  Time: %04d-%02d-%02d %02d:%02d:%02d%s", prefix, FRES_TIME->tm_year + 1900, FRES_TIME->tm_mon + 1, FRES_TIME->tm_mday, FRES_TIME->tm_hour, FRES_TIME->tm_min, FRES_TIME->tm_sec, g_linebreak_val);
	WriteParameters(FRES, spec, (char*)prefix);
}

/* Close result file and write parameters */
TARGET void CloseResFile()
{
	fclose(FRES);
}

/* Open temp files */
TARGET void OpenTempFiles(uint n, const char* spec, byte flag[])
{
	TEMP_NAMES = new char*[n];
	TEMP_FILES = new FILE*[n];
	TEMP_BUFS = new char*[n];

	for (uint i = 0; i < n; ++i)
	{
		if (flag && !flag[i]) continue;
		TEMP_BUFS[i] = new char[LINE_BUFFER];
		TEMP_NAMES[i] = new char[PATH_LEN];
		TEMP_FILES[i] = FOpen(TEMP_NAMES[i], "wb+", "%svcfpop_%04d_%02d_%02d_%02d_%02d_%02d%s%d%s",
			g_tmpdir_val.c_str(), FRES_TIME->tm_year + 1900, FRES_TIME->tm_mon + 1, FRES_TIME->tm_mday, FRES_TIME->tm_hour, FRES_TIME->tm_min, FRES_TIME->tm_sec, spec, i + 1, ".tmp");
		setvbuf(TEMP_FILES[i], TEMP_BUFS[i], _IOFBF, LINE_BUFFER);
	}
}

/* Close and Remove temp files */
TARGET void RemoveTempFiles(uint n)
{
	for (uint i = 0; i < n; ++i)
	{
		fclose(TEMP_FILES[i]);
		remove(TEMP_NAMES[i]);
		delete[] TEMP_NAMES[i];
		delete[] TEMP_BUFS[i];
	}
	delete[] TEMP_BUFS;
	delete[] TEMP_NAMES;
	delete[] TEMP_FILES;
}

/* Merge and close temp files */
TARGET void JoinTempFiles(uint n, byte flag[])
{
	int64 buflen = LINE_BUFFER;
	char* buf = new char[buflen];
	for (uint i = 0; i < n; ++i)
	{
		if (flag && !flag[i]) continue;
		fseeko64(TEMP_FILES[i], 0ll, SEEK_SET);
		while (!feof(TEMP_FILES[i]))
		{
			int64 t = fread(buf, 1, buflen, TEMP_FILES[i]);
			fwrite(buf, 1, t, FRES);
		}
		fclose(TEMP_FILES[i]);
		remove(TEMP_NAMES[i]);
		delete[] TEMP_NAMES[i];
		delete[] TEMP_BUFS[i];
	}
	delete[] TEMP_BUFS;
	delete[] TEMP_NAMES;
	delete[] TEMP_FILES;
	delete[] buf;
}

/* Move a file */
void FileMove(const char* src, const char* dst)
{
	rename(src, dst);
}

/* Does file exists */
bool FileExists(const char* file)
{
	return exists(file);
}

/* Delete a file */
void FileDelete(const char* file)
{
	remove(file);
}

/* Create a directory */
void DirectoryCreate(const char* path)
{
	create_directory(path);
}