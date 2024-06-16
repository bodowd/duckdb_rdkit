var duckdb = require('../../duckdb/tools/nodejs');
var assert = require('assert');

describe(`duckdb_rdkit extension`, () => {
    let db;
    let conn;
    before((done) => {
        db = new duckdb.Database(':memory:', {"allow_unsigned_extensions":"true"});
        conn = new duckdb.Connection(db);
        conn.exec(`LOAD '${process.env.DUCKDB_RDKIT_EXTENSION_BINARY_PATH}';`, function (err) {
            if (err) throw err;
            done();
        });
    });

    it('duckdb_rdkit function should return expected string', function (done) {
        db.all("SELECT duckdb_rdkit('Sam') as value;", function (err, res) {
            if (err) throw err;
            assert.deepEqual(res, [{value: "DuckdbRdkit Sam 🐥"}]);
            done();
        });
    });

    it('duckdb_rdkit_openssl_version function should return expected string', function (done) {
        db.all("SELECT duckdb_rdkit_openssl_version('Michael') as value;", function (err, res) {
            if (err) throw err;
            assert(res[0].value.startsWith('DuckdbRdkit Michael, my linked OpenSSL version is OpenSSL'));
            done();
        });
    });
});