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

#include "sqlite_wrapper.hpp"

namespace utilities
{
    //!
    //! Default constructor
    //!
    SqliteVariable::SqliteVariable()
        :
        m_name(""),
        m_type(""),
        m_data("0")
    {}    
    //!
    //! Copy constructor
    //!
    SqliteVariable::SqliteVariable(const SqliteVariable& other)
        :
        m_name(other.m_name),
        m_type(other.m_type),
        m_data(other.m_data)
    {}
    //!
    //! SqliteVariable Constructor without setting data values
    //!
    SqliteVariable::SqliteVariable(
    const std::string name,     //!<    SQL variable name
    const sql::sqlType_t type)       //!<    SQL variable type
        :
        m_name(name),
        m_data("0")
    {
        switch(type)
        {
            case sql::_TEXT_: 
                m_type = "TEXT";
                break;
            case sql::_INT_:
                m_type = "INT";
                break;
            case sql::_REAL_:
                m_type = "REAL";
                break;
            case sql::_BLOB_:
                m_type = "BLOB";
                break;
        }
    }
    //!
    //! Overload the comparison operator to declare objects equal if
    //! they share the same name and type
    //!
    bool SqliteVariable::operator==(const SqliteVariable& rhs) const
    {
        return (m_name == rhs.m_name && m_type == rhs.m_type);
    }
    //!
    //! Overload the comparison operator to declare objects equal if
    //! they share the same name and type
    //!
    bool SqliteVariable::operator==(SqliteVariable& rhs) const
    {
        return (m_name == rhs.m_name && m_type == rhs.m_type);
    }
    //!
    //! Overload the assignment operator to update only the class data
    //!
    SqliteVariable& SqliteVariable::operator=(const SqliteVariable& other)
    {
        m_data = other.m_data;
        return *this;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief Return a string that can be used to label the variable in a SQL
    //! table creation script
    //!
    //! \return SQL code segment for table creation of the stored variable
    //////////////////////////////////////////////////////////////////////////////////
    std::string SqliteVariable::TableCreationId() const
    {
        std::stringstream tmp;
        tmp.str("");
        tmp<<m_name<<" "<<m_type<<" ";
        return tmp.str();
    }
    //!
    //! Get the name parameter
    //!
    std::string SqliteVariable::GetName() const
    {
        return m_name;
    }
    //!
    //! Get the type parameter
    //!
    std::string SqliteVariable::GetType() const
    {
        return m_type;
    }
    //!
    //! \brief Get the data parameter
    //!
    std::string SqliteVariable::GetData() const
    {
        return m_data;
    }
    //!
    //! \brief Set the data parameter
    //!
    void SqliteVariable::SetData(const std::string newData)
    {
         m_data = newData;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function converts a set of data into a series of characters
    //!
    //!	The characters can then be passed to SQL to be interpreted as binary (BLOB)
    //!	data.
    //!
    //////////////////////////////////////////////////////////////////////////////////
    std::string SqliteVariable::Serialize(
        dcmplx* input,     //!<    Unserialized input array
        const int dim)     //!<    Size of input array
        const
    {
        std::stringstream output;
        output.str("");
        for(int i=0; i<dim; ++i)
        {
            output<<std::real(input[i])<<" "<<std::imag(input[i])<<" ";
        }
        return output.str().c_str();
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function converts a set of serialized data back into a regular 
    //! array
    //!
    //!	This function would be applied to descramble BLOB data written using the
    //!	serialize function
    //!
    //////////////////////////////////////////////////////////////////////////////////
    void SqliteVariable::Unserialize(
        char* serializedData,      //!<    Buffer containing serialized data
        dcmplx* output,            //!<    Unserialized output array
        const int dim)             //!<    Size of output array
        const
    {
        std::stringstream ss(serializedData);
        std::string buffer;
        std::vector<std::string> arrayElement;
        while (ss >> buffer)
        {
            arrayElement.push_back(buffer);
        }
        dcmplx *p_out;
        p_out=output;
        for(int i=0; i<2*dim; i+=2)
        {
            *(p_out)=dcmplx(atof(arrayElement[i].c_str()),atof(arrayElement[i+1].c_str()));
            ++p_out;
        }
        return;
    }

    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //!
    //! This function adds a new SQL data field to the list of fields that
    //! we can include in an SQL table
    //!
    void SqliteRow::AddField(
        const SqliteVariable& member) //!<    New SQL data field
    {
         m_fields.push_back(member);
    }
    //!
    //! This function updates an existing SQL data field 
    //!
    void SqliteRow::UpdateValue(
        const SqliteVariable& member) //!<    New SQL data field
    {
        bool found = false;
        for(unsigned int i=0; i<m_fields.size(); ++i)
        {
            if(member == m_fields[i])
            {
                m_fields[i] = member;
                found = true;
                break;
            }
        }
        if(!found)
        {
            std::cerr<<"\n\tSQLITE UPDATE MEMBER WARNING: existing field not found - did not update!"<<std::endl;   
        }
        return;
    }
    //!
    //! This function removes a SQL data field to the list of fields that
    //! we can include in an SQL table
    //!
    void SqliteRow::RemoveField(
        const SqliteVariable& member) //!<    SQL data field to remove
    {
        bool found = false;
        for(unsigned int i=0; i<m_fields.size(); ++i)
        {
            if(member == m_fields[i])
            {
                m_fields.erase(m_fields.begin()+i);
                found = true;
                break;
            }
        }
        if(!found)
        {
            std::cerr<<"\n\tSQLITE REMOVE MEMBER WARNING: existing field not found - did not remove!"<<std::endl;   
        }
        return;
    }
    //!
    //! \brief Assignment operator overload
    //!
    SqliteRow& SqliteRow::operator=(SqliteRow& lhs)
    {
        m_fields = lhs.m_fields;
        return *this;
    }
    //!
    //! \brief Assignment operator overload
    //!
    SqliteRow& SqliteRow::operator=(const SqliteRow& lhs)
    {
        m_fields = lhs.m_fields;
        return *this;
    }
    //!
    //! \brief Return the number of fields included
    //!
    unsigned int SqliteRow::GetLength() const
    {
        return m_fields.size();
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Return the list of fields currently associated with the class
    //!
    //! \return List of fields vector
    //!
    //////////////////////////////////////////////////////////////////////////////////   
    std::vector<SqliteVariable> SqliteRow::GetCurrentFields() const
    {
        return m_fields;
    }
    //!
    //! Public interface for the AddField function, passes variables
    //! to an SqliteVariable constructor
    //!
    void SqliteRow::AddField(const std::string name, const sql::sqlType_t type)
    {
        this->AddField(SqliteVariable(name, type));
    }
    //!
    //! Public interface for the UpdateValue function, passes variables
    //! to an SqliteVariable constructor
    //!
    void SqliteRow::UpdateValue(const std::string name,const sql::sqlType_t type)
    {
        this->UpdateValue(SqliteVariable(name,type));
    }
    
    //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //!
    //! Default constructor for the Sqlite class
    //!
    Sqlite::Sqlite()
    :
        m_db(0),
        m_dbOpen(sql::_NOT_OPEN_),
        m_timeUpdatedTrigger(false)
    {
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Constructor for the sqlite class
    //!
    //!	Opens a database with the file name databaseName
    //!
    //!	If the file already exists then it won't be overwritten, only appended
    //!
    //!	The db_open flag is then set to 1
    //////////////////////////////////////////////////////////////////////////////////
    Sqlite::Sqlite(
        const std::string fileName,      //!<   Name of Sqlite database
        const sql::sqlCreate_t flag)      //!<   Flag specifying whether to create
                                         //!    a new database or open an existing one   
    :
        m_db(0),
        m_dbOpen(sql::_NOT_OPEN_),
        m_timeUpdatedTrigger(false)
    {
        utilities::cout.MainOutput()<<"\n\tUsing SQLite version: "<<sqlite3_libversion()<<std::endl<<std::endl;
        std::stringstream dbName;
        dbName.str("");
        dbName<<fileName<<".sql";
        int errorFlag = 0;
        if(sql::_READ_EXISTING_ == flag)
        {
            errorFlag = sqlite3_open_v2(dbName.str().c_str(), &m_db, SQLITE_OPEN_READONLY, NULL);
        }
        else if(sql::_CREATE_NEW_ == flag)
        {
            errorFlag = sqlite3_open(dbName.str().c_str(),&m_db);
        }
        else if(sql::_EDIT_EXISTING_ == flag)
        {
            errorFlag = sqlite3_open_v2(dbName.str().c_str(), &m_db, SQLITE_OPEN_READWRITE, NULL);
        }
        if(errorFlag >= 1)
        {
            std::cerr<<"\n\tERROR WITH DATABASE: "<<dbName.str()<<std::endl;
            std::cerr<<"\n\tSQLITE MESSAGE:\n\n\t"<<sqlite3_errmsg(m_db)<<std::endl;
            sqlite3_close(m_db);
        }
        else
        {
            m_dbOpen = sql::_OPEN_;
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Destructor for the sqlite class
    //!
    //! Close the database if it's still open
    //////////////////////////////////////////////////////////////////////////////////
    Sqlite::~Sqlite()
    {
        if(sql::_OPEN_ == m_dbOpen) sqlite3_close(m_db);
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	A function to check if the database is successfully opened
    //!
    //! \return true if it is open, false if not
    //////////////////////////////////////////////////////////////////////////////////
    bool Sqlite::IsOpen() const
    {
        return (sql::_OPEN_ == m_dbOpen);
    }
    //!
    //! This is a wrapper function that passes a SQL script to the sqlite
    //! interface. Error messages are passed back to be printed out
    //!
    int Sqlite::ExecuteScript(
        const std::string sqlScript,    //!<    Sqlite script to be executed
        const bool getOutput)           //!<    Flag to retrieve output
    {
        if(sql::_OPEN_ == m_dbOpen)
        {
            char *errMsg;
            char** returnVal;
            int nRow;
            int nCol;
            int flag = sqlite3_get_table(m_db, sqlScript.c_str(), &returnVal, &nRow, &nCol, &errMsg);
            if(getOutput)
            {
                m_nRow = nRow;
                m_nCol = nCol;
                m_columnHead.clear();
                m_output.clear();
                for(int i=0; i < nCol; ++i)
                {
                   m_columnHead.push_back(returnVal[i]);   // First row heading //
                }
                for(int i=0; i < nCol*nRow; ++i)
                {
                   m_output.push_back(returnVal[nCol+i]);
                }
            }
            sqlite3_free_table(returnVal);
            if(0 != errMsg)
            {
                std::cerr<<"\n\n\t"<<sqlScript<<std::endl;
                std::cerr<<"\n\tSQLITE ERROR MESSAGE:\n\n\t\t"<<errMsg<<std::endl;
            }
            return flag;
        }
        else
        {
            std::cerr<<"\n\tSQLITE EXECUTE SCRIPT ERROR: Database not open!"<<std::endl;
        }
        return 1;
    }
    //!
    //! This function sets flags to include a time added trigger field 
    //! to the SQL table
    //!
    void Sqlite::AddTimeUpdatedTrigger(
        const std::string name)     //!<    SQL name of the date field
    {
        m_timeUpdatedTrigger = true;
        m_timeFieldName = name;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function generates and executes an SQLite script to create a
    //! new table in the opened database. An existing table of the same name is
    //! not overwritten
    //!
    //! \return true if there was an error, false otherwise
    //////////////////////////////////////////////////////////////////////////////////
    bool Sqlite::CreateTable(
        const std::string tableName,    //!<    Name of the table
        SqliteRow* fields)              //!<    List of fields to create
    {
        //	Set up the script to generate a new table from scratch, if required
        //  otherwise return the highest ID value from the existing table
        //  to allow appending
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<
        "CREATE TABLE IF NOT EXISTS "<<tableName<<"(	"
        "			   	Id INTEGER PRIMARY KEY,			";
        for(unsigned int i=0; i<fields->m_fields.size(); ++i)
        {
            sqlGenScript<<fields->m_fields[i].TableCreationId();
            if(i+1 < fields->m_fields.size())
            {
                sqlGenScript<<" , ";
            }
        }
        if(m_timeUpdatedTrigger)
        {
            sqlGenScript<<" , "<<m_timeFieldName<<" TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP";
        }
        sqlGenScript<<");";
        if(this->ExecuteScript(sqlGenScript.str(),false))
        {
            std::cerr<<"\n\tERROR WITH SQL CREATE TABLE"<<std::endl;
            return true;
        }
        if(m_timeUpdatedTrigger)
        {
            sqlGenScript.str("");
            sqlGenScript<<
            "CREATE TRIGGER IF NOT EXISTS update_"<<tableName<<"_trigger	     "
            "AFTER UPDATE ON "<<tableName<<" FOR EACH ROW						 "
            "BEGIN																 "
            "																	 "
            "UPDATE "<<tableName<<" SET "<<m_timeFieldName<<" = CURRENT_TIMESTAMP "
            "WHERE rowid = old.rowid;											 "
            "																	 "
            "END;																 "
            ;
            if(this->ExecuteScript(sqlGenScript.str(),false))
            {
                std::cerr<<"\n\tERROR WITH SQL CREATE TABLE"<<std::endl;
                return 0;
            }
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function generates and executes an SQlite script to remove an 
    //! entire Sql table
    //!
    //! \return true if there was an error, false otherwise
    //////////////////////////////////////////////////////////////////////////////////

    bool Sqlite::DeleteTable(
        const std::string tableName)        //!<    Name of the table
    {
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<
        "DROP TABLE "<<tableName;
        if(this->ExecuteScript(sqlGenScript.str(),false))
        {
            std::cerr<<"\n\tERROR WITH SQL DELETE TABLE"<<std::endl;
            return true;
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function generates and executes an Sqlite script to insert a
    //! new set of entries into an existing table
    //!
    //! \return true if there was an error, false otherwise
    //////////////////////////////////////////////////////////////////////////////////
    bool Sqlite::InsertIntoTable(
        const std::string tableName,       //!<    Name of the table
        std::vector<SqliteRow>* data)      //!<  List of data to be entered
    {
        const int nbrRows = data->size();
        const int nbrFields = (*data)[0].m_fields.size();
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<std::setprecision(stringStreamPrecision)<<"INSERT INTO "<<tableName<<"(";
        for(int i=0; i<nbrFields; ++i)
        {
            sqlGenScript<<(*data)[0].m_fields[i].GetName();
            if(i+1 < nbrFields)
            {
                sqlGenScript<<" , "; 
            }
        }
        sqlGenScript<<") VALUES ";
        for(int j=0; j<nbrRows; ++j)
        {	
            sqlGenScript<<"(";
            SqliteRow row = (*data)[j];
            for(int i=0; i<nbrFields; ++i)
            {
                if("TEXT" == row.m_fields[i].GetType() || "BLOB" == row.m_fields[i].GetType())
                {
                    //  Put in enclosing '' for TEXT or BLOB types
                    sqlGenScript<<"'"<<row.m_fields[i].GetData()<<"'";
                }
                else
                {
                    sqlGenScript<<row.m_fields[i].GetData();
                }
                if(i+1 < nbrFields)
                {
                    sqlGenScript<<" , "; 
                }
            }
            sqlGenScript<<")";
            if(j+1 < nbrRows)
            {
                sqlGenScript<<" , "; 
            }
        }
        if(this->ExecuteScript(sqlGenScript.str(),false))
        {
            std::cerr<<"\n\tERROR WITH SQL INSERT INTO TABLE"<<std::endl;
            return true;
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function generates and executes an Sqlite script to insert a
    //! single row into an sqlite table (this method is very inefficient - it's much
    //! better to insert all required rows in one command using a std::vector of
    //! SqliteRow objects
    //!
    //! \return true if there was an error, false otherwise
    //////////////////////////////////////////////////////////////////////////////////
    bool Sqlite::InsertIntoTable(const std::string tableName, SqliteRow* data)
    {
        std::vector<SqliteRow> temp(1);
        temp[0] = *data;
        return InsertIntoTable(tableName, &temp);
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function generates and executes an SQlite script to remove an
    //! existing entry into an existing table
    //!
    //! \return true if there was an error, false otherwise
    //////////////////////////////////////////////////////////////////////////////////         
    bool Sqlite::RemoveTableEntry(
        const std::string tableName,        //!<    Name of the table
        SqliteRow* fields)                  //!<    List of fields
    {
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<std::setprecision(stringStreamPrecision)<<
        "DELETE FROM "<<tableName<<" WHERE ";
        for(unsigned int i=0; i<fields->m_fields.size(); ++i)
        {
            sqlGenScript<<fields->m_fields[i].GetName()<<"="<<fields->m_fields[i].GetData();
            if(i+1 < fields->m_fields.size())
            {
                sqlGenScript<<" AND "; 
            }
        }
        if(this->ExecuteScript(sqlGenScript.str().c_str(),false))
        {
            std::cerr<<"\n\tERROR WITH SQL REMOVE TABLE ENTRY "<<std::endl;
            return true;
        }
        return false;
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Update an existing table entry
    //!
    //! \return true if there was an error, false otherwise
    //////////////////////////////////////////////////////////////////////////////////
    bool Sqlite::UpdateTableEntry(
        const std::string tableName,        //!<    Name of table
        const int sqlId,                    //!<    Specify the table row to update
        const SqliteVariable newValue)      //!<    Field to update and new value
    {
        if(sqlId<=0)
        {
            std::cerr<<"ERROR in UpdateTableEntry: sql ID <0 "<<std::endl;
            return true;
        }
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<std::setprecision(stringStreamPrecision)<<"UPDATE "<<tableName
        <<" SET "<<newValue.GetName()<<"="<<newValue.GetData()<<" WHERE id="<<sqlId<<";";
        if(this->ExecuteScript(sqlGenScript.str(),false))
        {
            std::cerr<<"\n\tERROR WITH SQL UPDATE TABLE ENTRY"<<std::endl;
            return true;
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Return selected values from the table matching a particular pattern
    //!
    //! \return A vector of fields list
    //////////////////////////////////////////////////////////////////////////////////
    std::vector<SqliteRow> Sqlite::RetrieveFromTable(
        const std::string tableName,//!<    Name of table
        SqliteRow* fields)          //!<    Search fields to match
    {
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<std::setprecision(stringStreamPrecision)<<
        "SELECT * FROM "<<tableName<<" WHERE ";
        for(unsigned int i=0; i<fields->m_fields.size(); ++i)
        {
            sqlGenScript<<fields->m_fields[i].GetName()<<"="<<fields->m_fields[i].GetData();
            if(i+1 < fields->m_fields.size())
            {
                sqlGenScript<<" AND "; 
            }
        }
        if(this->ExecuteScript(sqlGenScript.str().c_str(), true))
        {
            std::cerr<<"\n\tERROR WITH SQL RETRIEVE FROM TABLE"<<std::endl;
        }
        std::vector<SqliteRow> tmp(m_nRow);
        for(int i=0; i < m_nRow; ++i)
        {
            SqliteRow tmp1 = *fields;
            for(unsigned int j=1; j <= tmp1.m_fields.size(); ++j)
            {
                utilities::cout.MainOutput()<<m_output[i*m_nCol+j]<<std::endl;
                tmp1.m_fields[j-1].SetData(m_output[i*m_nCol+j]);
            }
            tmp[i] = tmp1;
        }
        return tmp;
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Return selected values from the table for a given SQL id 
    //! and put them in the m_fields list
    //!
    //! \return Sqlite table row data
    //////////////////////////////////////////////////////////////////////////////////
         
    SqliteRow Sqlite::RetrieveIdFromTable(
        const std::string tableName,           //!<    Name of table
        const int sqlId)                       //!<    Unique Id of the desired table entry
    {     
        SqliteRow row;
        if(sqlId<=0)
        {
            std::cerr<<"ERROR in UpdateTableEntry: sql ID <0 "<<std::endl;
            return row;
        }
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<std::setprecision(stringStreamPrecision)<<
        "SELECT * FROM "<<tableName<<" WHERE Id="<<sqlId;
        if(this->ExecuteScript(sqlGenScript.str().c_str(), true))
        {
            std::cerr<<"\n\tERROR WITH SQL RETRIEVE FROM TABLE"<<std::endl;
            return row;
        }
        for(int j=1; j < m_nCol; ++j)
        {
            SqliteVariable temp;
            temp.SetValues(m_columnHead[j], m_output[j]);
            row.m_fields.push_back(temp);
        }
        return row;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Print out the currently stored SQL table
    //!
    //! \return true if there was an error, false otherwise
    //////////////////////////////////////////////////////////////////////////////////    
    bool Sqlite::PrintTable(
        const std::string tableName)        //!<    Name of the table
    {
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<"SELECT * FROM "<<tableName;
        if(this->ExecuteScript(sqlGenScript.str().c_str(), true))
        {
            std::cerr<<"\n\tERROR WITH SQL PRINT TABLE"<<std::endl;
            return true;
        }
        int widths[m_nCol];
        for(int j=0; j < m_nCol; ++j)
        {
           widths[j] = std::max(m_output[j].length(),m_columnHead[j].length()) + 4; //  +4 for tab spacing
        }
        for(int j=0; j < m_nCol; ++j)
        {
           utilities::cout.MainOutput()<<std::setw(widths[j])<<std::left<<m_columnHead[j];
        }
        utilities::cout.MainOutput()<<std::endl;
        for(int i=0; i < m_nRow; ++i)
        {
            for(int j=0; j < m_nCol; ++j)
            {
                utilities::cout.MainOutput()<<std::setw(widths[j])<<std::left<<m_output[i*m_nCol+j];
            }
            utilities::cout.MainOutput()<<std::endl;
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Print out the currently stored SQL table to a text file
    //!
    //! \return true if there was an error, false otherwise
    //!
    //////////////////////////////////////////////////////////////////////////////////
    bool Sqlite::CopyTableToTextFile(
        const std::string tableName,        //!<    Name of the table
        const std::string fileName)         //!<    Name of text file

    {
        std::stringstream name;
        name.str("");
        name<<fileName<<".txt";
        std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
        std::ofstream   fout(name.str().c_str());
        std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
        bool flag = this->PrintTable(tableName);
        std::cout.rdbuf(cout_sbuf); // restore the original cout buffer
        return flag;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	Use a text file print out of an SQL table to generate a new sql
    //! table
    //!
    //! \return true if there was an error, false otherwise
    ////////////////////////////////////////////////////////////////////////////////// 
    bool Sqlite::CopyTextFileToTable(
        const std::string tableName,        //!<    Name of table
        const std::string fileName)         //!<    Text file name
    {
        //  TODO complete this function

        #if 0
        std::ifstream f_in;
        f_in.open(fileName.c_str());
        if(!f_in.is_open())
        {
            std::cerr<<"\n\tERROR WITH SQL CopyTextFileToTable: file "<<fileName<<" not found."<<std::endl;
            return true;
        }
        //  Clear current field data held by class
        m_fields.clear();
        //  First line contains the table headings
        std::string line;
        std::string field;
        std::getline(f_in,line);
        std::stringstream headings(line);
        utilities::cout.MainOutput()<<"\n\t - RECONSTRUCTING SQL TABLE FROM TEXT FILE. "; 
        utilities::cout.MainOutput()<<"PLEASE DEFINE DATA TYPES ON REQUEST:\n\n"
        utilities::cout.MainOutput()<<"\tOPTIONS: 0 for TEXT, 1 for INT, 2 for REAL, 3 for BLOB\n"<<std::endl;
        while(std::getline(headings,field))
        {
            if(field != "Id")
            {
                utilities::cout.MainOutput()<<"\t FIELD NAME: "<<field<<std::endl;
                int choice;
                std::cin>>choice;
                switch(choice)
                {
                    case 0:
                        this->AddField(field,sql::_TEXT_);
                        break;
                    case 1:
                        this->AddField(field,sql::_INT_);
                        break;
                    case 2:
                        this->AddField(field,sql::_REAL_);
                        break;
                    case 3:
                        this->AddField(field,sql::_BLOB_);
                        break;
                }
            
            }
        }
        std::vector<std::vector<SqliteVariable> > data;
        while(std::getline(headings,line))
        {}
        f_in.close();
        #endif
        return false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief  Get the maximum ID used in the current table
    //!
    //! \return Maximum table ID, or return 0 if the table is empty
    ////////////////////////////////////////////////////////////////////////////////// 
    int Sqlite::GetMaxId(
        const std::string tableName)        //!<    Name of table
    {
        std::stringstream sqlGenScript;
        sqlGenScript.str("");
        sqlGenScript<<std::setprecision(stringStreamPrecision)<<
        "SELECT COUNT(*) FROM "<<tableName;
        if(this->ExecuteScript(sqlGenScript.str().c_str(),true))
        {
            std::cerr<<"\n\tERROR WITH SQL GetMaxId"<<std::endl;
        }
        if(std::atoi(m_output[0].c_str())>0)
        {
            sqlGenScript.str("");
            sqlGenScript<<std::setprecision(stringStreamPrecision)<<
            "SELECT MAX(id) AS max_id FROM "<<tableName;
            if(this->ExecuteScript(sqlGenScript.str().c_str(),true))
            {
                std::cerr<<"\n\tERROR WITH SQL GetMaxId"<<std::endl;
            }
            return std::atoi(m_output[0].c_str());
        }
        else
        {
            return 0;
        }
    }
}   //  End namespace utilities
