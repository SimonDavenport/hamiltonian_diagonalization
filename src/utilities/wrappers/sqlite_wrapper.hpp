////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!		This file contains a bunch of functions used for interfacing with
//!		an SQLite database. Program source:
//!		http://souptonuts.sourceforge.net/code/simplesqlite3cpp2.cc.html
//!
//!                    Copyright (C) Simon C Davenport
//!                                                                             
//!     This program is free software: you can redistribute it and/or modify
//!     it under the terms of the GNU General Public License as published by
//!     the Free Software Foundation, either version 3 of the License,
//!     or (at your option) any later version.
//!                                                                             
//!     This program is distributed in the hope that it will be useful, but
//!     WITHOUT ANY WARRANTY; without even the implied warranty of
//!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//!     General Public License for more details.
//!                                                                             
//!     You should have received a copy of the GNU General Public License
//!     along with this program. If not, see <http://www.gnu.org/licenses/>.
//!                                                                             
////////////////////////////////////////////////////////////////////////////////

#ifndef _SQLITE_WRAPPER_HPP_INCLUDED_
#define _SQLITE_WRAPPER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../general/dcmplx_type_def.hpp"
#include "../general/template_tools.hpp"
#include "../general/cout_tools.hpp"
#include <sqlite3.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cassert> 
#if _DEBUG_
#include "../general/debug.hpp"
#endif
//////      ENUM DECLARATIONS       ////////////////////////////////////////
enum sqlFlag_t {_OPEN_=0,_NOT_OPEN_=1};            //!<   SQL Database open flag
enum sqlType_t {_TEXT_,_INT_,_REAL_,_BLOB_};       //!<   Possible SQL data types
enum sqlCreate_t {_CREATE_NEW_,_READ_EXISTING_,_EDIT_EXISTING_};   
                                    //!<   Flag to pass to the Sqlite constructor
namespace utilities
{
    class SqliteVariable
    {
        private:
        std::string m_name;         //!<    A string containing the name of the
                                    //!     SQL variable
        std::string m_type;         //!<    A string containing the type of the
                                    //!     SQL variable
        std::string m_data;         //!<    Data associated with this member
        public:
        SqliteVariable();
        SqliteVariable(const SqliteVariable& other);
        bool operator==(const SqliteVariable& rhs) const;
        bool operator==(SqliteVariable& rhs) const;
        SqliteVariable& operator=(const SqliteVariable& other);
        SqliteVariable(const std::string name, const sqlType_t type);

        //!
        //! Template method setting explicit data values
        //!
        template<typename T> 
        void SetValues(const std::string name, const T data)
        {
            m_name = name;
            //  Convert template typename to a string containing 
            //  the SQL type name and convert the input data into a string
            std::stringstream tmp;
            tmp.str("");
            tmp<<data;   
            m_data = tmp.str();
            if(is_same<T, int>::value)
            {     
                //  Set INT type if a type is not already defined by the constructor
                if("" == m_type) m_type = "INT";
                assert(m_type =="INT");
            }
            else if(is_same<T, double>::value)
            {
                //  Set REAL type if a type is not already defined by the constructor
                if("" == m_type) m_type = "REAL";
                assert(m_type =="REAL");
            }
            else
            {
                //  Set TEXT type if a type is not already defined by the constructor
                if("" == m_type) m_type = "TEXT";
                //  Allow BLOB type also here
                assert(m_type =="TEXT" || m_type == "BLOB");
            }
        }
        std::string TableCreationId() const;
        std::string GetName() const;
        std::string GetType() const;
        std::string GetData() const;
        void SetData(const std::string newData);
        std::string Serialize(dcmplx* input, const int dim) const;
        void Unserialize(char* serializedData, dcmplx* output, const int dim ) const;
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A class to contain is list of SqliteVariables to form a single
    //! row in an Sqlite table
    ////////////////////////////////////////////////////////////////////////////////
    
    class SqliteRow
    {
        private:
        std::vector<SqliteVariable> m_fields;       //!<    List of SQL data names and types
        void AddField(const SqliteVariable& member);
        void UpdateValue(const SqliteVariable& member);
        void RemoveField(const SqliteVariable& member);
        friend class Sqlite;
        public:
        SqliteRow& operator=(SqliteRow& lhs);
        SqliteRow& operator=(const SqliteRow& lhs);
        unsigned int GetLength() const;
        std::vector<SqliteVariable> GetCurrentFields() const;
        void AddField(const std::string name,const sqlType_t type);
        //!
        //! Template for the AddField function to pass arguments to the 
        //! SqliteVariable constructor
        //!
        template<typename T>
        void AddField(const std::string name, const T data)
        {
            SqliteVariable temp;
            temp.SetValues(name,data);
            this->AddField(temp);
        }
        //!
        //! Template for the UpdateValue function to pass arguments to the 
        //! SqliteVariable constructor
        //!
        template<typename T>
        void UpdateValue(const std::string name, const sqlType_t type, const T data)
        {
            SqliteVariable temp(name,type);
            temp.SetValues(name,data);
            this->UpdateValue(temp);
        }
        //!
        //! Template for the UpdateValue function to pass arguments to the 
        //! SqliteVariable constructor
        //!
        template<typename T>
        void UpdateValue(const std::string name, const T data)
        {
            SqliteVariable temp;
            temp.SetValues(name, data);
            this->UpdateValue(temp);
        }
        
        void UpdateValue(const std::string name, const sqlType_t type);
        //!
        //! Template to get the parameter value corresponding to the specified 
        //! SQL field and type from the list of fields sorted in m_data
        //!
        template<typename T> 
        bool GetValue(
            const std::string name,     //!<    Name of SQL field
            T& value)                   //!<    Data value to be returned
        {
            for(int i=0; i<m_fields.size(); ++i)
            {
                if(name == m_fields[i].GetName())
                {
                    std::stringstream convert;
                    convert.str("");
                    convert<<(m_fields[i].GetData());
                    convert >> value;
                    return false;
                }
            }
            std::cout<<"\n\tSQLITE GetValue WARNING: existing field "<<name<<" not found!"<<std::endl; 
            return true;
        }
    };

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A wrapper class for SQLlite functions
    //!
    ////////////////////////////////////////////////////////////////////////////////
    class Sqlite
    {
	    private:
        const static int stringStreamPrecision=15;  //!<    Precision when writing
                                                    //!     numbers in a string stream
        sqlite3* m_db;              //!<    Pointer to SQL database object
        sqlFlag_t m_dbOpen;         //!<    Flag to specify if the database is opened
        bool m_timeUpdatedTrigger;  //!<    Flag to specify a time field
        std::string m_timeFieldName;//!<    Name of SQL time field
        std::vector<std::string> m_columnHead;
        std::vector<std::string> m_output;
        int m_nRow;               //!<    number of rows in output
        int m_nCol;               //!<    number of columns in output
        int ExecuteScript(const std::string script, const bool getOutput);

	    public:
        Sqlite();
        Sqlite(const std::string fileName, const sqlCreate_t flag);
        ~Sqlite();
        bool IsOpen() const;
        void AddTimeUpdatedTrigger(const std::string name);
        bool CreateTable(const std::string tableName, SqliteRow* fields);
        bool DeleteTable(const std::string tableName);
        int GetMaxId(const std::string tableName);
        bool InsertIntoTable(const std::string tableName, std::vector<SqliteRow>* data);
        bool InsertIntoTable(const std::string tableName, SqliteRow* data);
        bool RemoveTableEntry(const std::string tableName, SqliteRow* fields);
        bool UpdateTableEntry(const std::string tableName, const int sqlId,
                              const SqliteVariable newValue);
        SqliteRow RetrieveIdFromTable(const std::string tableName, const int sqlId);
        std::vector<SqliteRow> RetrieveFromTable(const std::string tableName,
                                                 SqliteRow* fields);
        bool PrintTable(const std::string tableName);
        bool CopyTableToTextFile(const std::string tableName, const std::string fileName);
        bool CopyTextFileToTable(const std::string tableName, const std::string fileName);
    };
}   //  End namespace utilities
#endif
