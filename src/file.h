/* File IO */

#pragma once
#include "vcfpop.h"

/* Get file size */
TARGET int64 GetFileLen(char* file);

/* Read a T from file */
template<typename T>
TARGET T FGet(FILE* file, T buf)
{
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, sizeof(T));
	else
		fread(&buf, sizeof(T), file);
	return buf;
}

/* Read bytes from file */
template<typename T>
TARGET T* FGet(FILE* file, T* buf, uint64 size)
{
	if (g_format_val < 0)
		gzread((gzFile)file, buf, sizeof(T) * (uint)size);
	else
		fread(buf, 1, sizeof(T) * size, file);
	return buf;
}

/* Read a double from file */
TARGET double FGetDouble(FILE* file);

/* Read a float from file */
TARGET float FGetFloat(FILE* file);

/* Read an int64 from file */
TARGET int64 FGetInt64(FILE* file);

/* Read an uint64 from file */
TARGET uint64 FGetUint64(FILE* file);

/* Read an int from file */
TARGET int FGetInt(FILE* file);

/* Read an uint from file */
TARGET uint FGetUint(FILE* file);

/* Read a short from file */
TARGET short FGetShort(FILE* file);

/* Read an ushort from file */
TARGET ushort FGetUshort(FILE* file);

/* Read a byte from file */
TARGET byte FGetByte(FILE* file);

/* Read a char from file */
TARGET char FGetChar(FILE* file);

/* Read a typed int from bcf file */
TARGET uint FGetTypedInt(FILE* file);

/* Is end of file */
TARGET bool FEof(FILE* file);

/* Get compressed data offset */
TARGET int64 FOffset(FILE* file);

/* Get decompressed data offset */
TARGET int64 FTell(FILE* file);

/* Read decompressed data */
TARGET int64 FRead(void* buf, uint64 element, uint64 count, FILE* file);

/* Change file offset for decompressed data */
TARGET int64 FSeek(FILE* file, int64 offset, int origin);

/* Open an decompressed or compressed file */
TARGET FILE* FOpen(const char* file, const char* spec);

/* Close an decompressed or compressed file */
TARGET int FClose(FILE* file);

/* Format file name and open file */
TARGET FILE* FOpen(char* buf, const char* type, const char* fmt ...);

/* Open result file and write parameters */
TARGET void OpenResFile(const char* _spec, const char* title);

/* Close result file and write parameters */
TARGET void CloseResFile();

/* Open temp files */
TARGET void OpenTempFiles(uint n, const char* spec, byte flag[] = NULL);

/* Merge and close temp files */
TARGET void JoinTempFiles(uint n, byte flag[] = NULL);

/* Move a file */
void FileMove(const char* src, const char* dst);

/* Does file exists */
bool FileExists(const char* file);

/* Delete a file */
void FileDelete(const char* file);

/* Create a directory */
void DirectoryCreate(const char* path);