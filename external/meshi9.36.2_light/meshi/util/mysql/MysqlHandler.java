package meshi.util.mysql;
import java.io.*;
import java.util.*;
import java.sql.*;
import meshi.util.*;

/***
    This class accesses MSQL-based database in order to update the fragments table.
 ***/
public class MysqlHandler {
    private Connection conn = null;

    private static final int DEBUG_LEVEL = 1;
    private static final int RC_NO_SUCH_TABLE = 1146;

    /*** get username/password from the command file ***/
    public MysqlHandler(CommandList commands)
        {
           String user     = commands.firstWord("mysqluser").secondWord();
           String password = commands.firstWord("mysqlpassword").secondWord();
           String database = commands.firstWord("mysqldatabase").secondWord();
           
           try
           {
               // open connection
               Class.forName("com.mysql.jdbc.Driver").newInstance();
               conn = DriverManager.getConnection("jdbc:mysql:///"+database,
               user, password);

               if(!conn.isClosed())
                 System.out.println("Successfully connected to " +
                 "MySQL server using TCP/IP...");
             }
           catch(Exception e)
           {
               System.out.println("Failed to connect to " +
                 "MySQL server using TCP/IP...");
               conn = null;
               if (DEBUG_LEVEL>=1)
                   e.printStackTrace();
           }
        }

    /*** get username/password explicitly ***/
    public MysqlHandler(String user, String password, String database)
        {
           try
           {
               // open connection
               Class.forName("com.mysql.jdbc.Driver").newInstance();
               conn = DriverManager.getConnection("jdbc:mysql:///"+database,
               user, password);

               if(!conn.isClosed())
                 System.out.println("Successfully connected to " +
                 "MySQL server using TCP/IP...");
             }
           catch(Exception e)
           {
               e.printStackTrace();
           }
        }
    
    /***********************************************************************************/
    /***********************       VERIFICATION METHODS      ***************************/
    /***********************************************************************************/
    public boolean connected() { return (conn != null); }

    /***********************************************************************************/
    /***********************       TABLE CREATION METHODS    ***************************/
    /***********************************************************************************/    
    public boolean createTable(String tableName, String tableParameters)
        {
            /*** drop table if exists ***/
            try
            {
            Statement s = conn.createStatement ();
            String dropQuery = "DROP TABLE "+tableName+";";
            s.executeUpdate (dropQuery);
            ResultSet rs = s.getResultSet ();
            s.close();
            }
            catch (SQLException e)
            {
            }

            /*** create table anew ***/
            try
            {
              Statement s = conn.createStatement ();
              String createQuery = "CREATE TABLE "+tableName+tableParameters;
              s.executeUpdate (createQuery);
              ResultSet rs = s.getResultSet ();
              s.close();
            }
            catch (SQLException e)
            {
                e.printStackTrace();
                return false;
            }
            return true;
        }

    public boolean createNonExistingTable(String tableName, String tableParameters)
        {
            if (tableName == null)
                return false;
            
            /*** check if table exists ***/
            try
            {
               Statement s = conn.createStatement ();
               String dropQuery = "SELECT COUNT(*) FROM "+tableName+";";
               s.executeQuery (dropQuery);
               ResultSet rs = s.getResultSet ();
               s.close();
               System.out.println("Table "+tableName+" exists, there is no need to create it");
               return false;
            }
            catch (SQLException e)
            {
                /*** we could be here because the table does not exist ***/
                try
                 {
                    Statement s = conn.createStatement ();
                    String createQuery = "CREATE TABLE "+tableName+tableParameters;
                    s.executeUpdate (createQuery);
                    ResultSet rs = s.getResultSet ();
                    s.close();
                    return true;
                 }
                catch (SQLException e1)
                 {
                    e1.printStackTrace();
                    return false;
                  }
            }

        }

    /***********************************************************************************/
    /***********************        TABLE INSERT METHODS     ***************************/
    /***********************************************************************************/
    /***
        Insert pre-setResidue values
     ***/
    public boolean insertIntoTable(String tableName, String tableParameters, String values)
        {            
            /*** intsert into table  ***/
            try
            {
              Statement s = conn.createStatement ();
              String insertQuery = "INSERT INTO "+tableName+" "+tableParameters+" VALUES "+values+";";
              if (DEBUG_LEVEL >=2)
                  System.out.println("Query: "+insertQuery);
              s.executeUpdate (insertQuery);
              ResultSet rs = s.getResultSet ();
              s.close();
            }
            catch (SQLException e)
            {
                e.printStackTrace();
                return false;
            }
            return true;
        }

    /***********************************************************************************/
    /***********************          TABLE EXISTENCE        ***************************/
    /***********************************************************************************/
    public boolean tableExists(String tableName)
        {
            if (tableName != null)
            {
                /*** with a simple query verify table's existence ***/
                try
                {
                  // check for table FRAGMENTS
                  Statement s = conn.createStatement ();
                  s.executeQuery ("SELECT COUNT(*) FROM "+tableName);
                  ResultSet rs = s.getResultSet ();
                  s.close();
                  return true;
                }
                catch (SQLException e)
                {
                }
            }
            return false;
        }

    /***********************************************************************************/
    /***********************        TABLE SELECT METHODS     ***************************/
    /***********************************************************************************/
    /*** select all ***/
    public ArrayList selectAll(String tableName, int numOfAttributes)
        {
            try
            {
              Statement s = conn.createStatement ();
              String selectQuery = "SELECT * FROM "+tableName;
              
              if (DEBUG_LEVEL >=3)
                  System.out.println("Query: "+selectQuery);
              
              s.executeQuery (selectQuery);
              ResultSet rs = s.getResultSet ();
              ArrayList list = new ArrayList();
              
              while (rs.next())
                {
                    String row[] = new String[numOfAttributes];
                    for (int i=0; i<numOfAttributes; i++)
                      row[i] = rs.getString(i+1);
                    list.add(row);
                }
              
              s.close();              
              return list;
            }
            catch (SQLException e)
            {
                e.printStackTrace();
                return null;
            }
        }

    /*** select row count, with condition ***/
    public int selectCount(String tableName, String whereCondition)
        {
            try
            {
              Statement s = conn.createStatement ();
              String selectQuery = "SELECT COUNT(*) as cnt FROM "+tableName+" WHERE "+whereCondition;
              if (DEBUG_LEVEL >=3)
                  System.out.println("Query: "+selectQuery);
              s.executeQuery (selectQuery);
              ResultSet rs = s.getResultSet ();
              
              int rc = -1;
              if (rs.next())
                rc = rs.getInt("cnt");
              s.close();
              return rc;
            }
            catch (SQLException e)
            {
                e.printStackTrace();
                return -1;
            }
        }

    /*** select one attribute, with condition ***/
    public double selectDouble(String tableName, String whereCondition, String selectAttribute)
        {
            try
            {
              Statement s = conn.createStatement ();
              String selectQuery = "SELECT "+selectAttribute+" FROM "+tableName+" WHERE "+whereCondition;
              if (DEBUG_LEVEL >=3)
                  System.out.println("Query: "+selectQuery);
              s.executeQuery (selectQuery);
              ResultSet rs = s.getResultSet ();

              double w = -1;
              if (rs.next())
               w = rs.getDouble(selectAttribute);
              s.close();
              return w;
            }
            catch (SQLException e)
            {
                e.printStackTrace();
                return -1;
            }
        }

    /*** select one attribute, with condition, return a list ***/
    public ArrayList selectDoubleList(String tableName, String whereCondition, String selectCondition, String selectAttribute)
        {
            try
            {
              Statement s = conn.createStatement ();
              String selectQuery = "SELECT "+selectCondition+" FROM "+tableName+" WHERE "+whereCondition;
              
              if (DEBUG_LEVEL >=3)
                  System.out.println("Query: "+selectQuery);
              
              s.executeQuery (selectQuery);
              ResultSet rs = s.getResultSet ();
              ArrayList list = new ArrayList();
              
              while (rs.next())
              {
                  list.add(new Double(rs.getDouble(selectAttribute)));
              }

              if (DEBUG_LEVEL >=3)
                  System.out.println("found "+list.size()+" doubles");
              
              s.close();              
              return list;
            }
            catch (SQLException e)
            {
                e.printStackTrace();
                return null;
            }
        }
    
    /***********************************************************************************/
    /***********************        GENERIC METHODS          ***************************/
    /***********************************************************************************/
    /*****
          Close connection to the MYSQL server
     *****/
    public void close()
        {
            try {
                 if(conn != null)
                    conn.close();
                 System.out.println("Successfully disconnected from MySQL server...");
                
               } catch(SQLException e) {}
            catch(Exception e)
             {
               System.out.println("Exception: " + e.getMessage());
               e.printStackTrace();
             } 
        }
}
