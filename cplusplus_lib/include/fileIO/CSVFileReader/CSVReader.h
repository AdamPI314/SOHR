#ifndef __CSVREADER_H_
#define __CSVREADER_H_
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility> //for std::pair
#include <iterator> //for std::vector<T>::iterator
#include <vector>
#include <algorithm> //for std:copy
#include <sstream>


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This class is written to deal with read/write file stuff.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//look the website below for more info
//http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c

// Best read CSV file, could read other files with small modification
//I would just create a class representing a row.
//Then stream into that object:
class CSVRow
{
public:
	static std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str);
public:
	std::string const& operator[](std::size_t index) const
	{
		return m_data[index];
	}

	std::size_t size() const
	{
		return m_data.size();
	}

	void readNextRow(std::istream& str);

private:
	std::vector<std::string>    m_data;
};

//But with a little work we could technically create an iterator:
class CSVIterator
{
public:
	typedef std::input_iterator_tag     iterator_category;
	typedef CSVRow                      value_type;
	typedef std::size_t                 difference_type;
	typedef CSVRow*                     pointer;
	typedef CSVRow&                     reference;

public:
	CSVIterator(std::istream& str) :m_str(str.good() ? &str : NULL) { ++(*this); }
	CSVIterator() :m_str(NULL) {}

	//// Pre Increment
	CSVIterator& operator++() { if (m_str) { m_row.readNextRow(*m_str); m_str = m_str->good() ? m_str : NULL; }return *this; }
	//// Post increment
	CSVIterator operator++(int) { CSVIterator tmp(*this); ++(*this); return tmp; }
	CSVRow const& operator*()   const { return m_row; }
	CSVRow const* operator->()  const { return &m_row; }

	bool operator==(CSVIterator const& rhs) { return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL))); }
	bool operator!=(CSVIterator const& rhs) { return !((*this) == rhs); }
private:
	std::istream*       m_str;
	CSVRow              m_row;
};


#endif
