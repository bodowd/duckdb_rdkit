import duckdb
import os
import pytest

# Get a fresh connection to DuckDB with the duckdb_rdkit extension binary loaded
@pytest.fixture
def duckdb_conn():
    extension_binary = os.getenv('DUCKDB_RDKIT_EXTENSION_BINARY_PATH')
    if (extension_binary == ''):
        raise Exception('Please make sure the `DUCKDB_RDKIT_EXTENSION_BINARY_PATH` is set to run the python tests')
    conn = duckdb.connect('', config={'allow_unsigned_extensions': 'true'})
    conn.execute(f"load '{extension_binary}'")
    return conn

def test_duckdb_rdkit(duckdb_conn):
    duckdb_conn.execute("SELECT duckdb_rdkit('Sam') as value;");
    res = duckdb_conn.fetchall()
    assert res[0][0] == "DuckdbRdkit Sam 🐥"

def test_duckdb_rdkit_openssl_version_test(duckdb_conn):
    duckdb_conn.execute("SELECT duckdb_rdkit_openssl_version('Michael');");
    res = duckdb_conn.fetchall()
    assert res[0][0][0:51] == "DuckdbRdkit Michael, my linked OpenSSL version is OpenSSL"
