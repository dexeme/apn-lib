PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS apn_family (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    description TEXT
);

CREATE TABLE IF NOT EXISTS apn_function (
    id INTEGER PRIMARY KEY,
    local_id INTEGER,
    dimension INTEGER NOT NULL,
    canonical_expression TEXT NOT NULL,
    canonical_json TEXT,
    normalized_support TEXT,
    equivalent_to TEXT,
    walsh_spectrum TEXT,
    source_label TEXT,
    notes TEXT,
    UNIQUE (dimension, local_id)
);

CREATE TABLE IF NOT EXISTS apn_function_representation (
    id INTEGER PRIMARY KEY,
    function_id INTEGER NOT NULL,
    representation_type TEXT NOT NULL,
    representation_json TEXT NOT NULL,
    display_expression TEXT,
    parameters TEXT,
    FOREIGN KEY (function_id) REFERENCES apn_function(id)
);

CREATE TABLE IF NOT EXISTS apn_function_family (
    function_id INTEGER NOT NULL,
    family_id INTEGER NOT NULL,
    relation TEXT NOT NULL,
    parameters TEXT,
    PRIMARY KEY (function_id, family_id, relation),
    FOREIGN KEY (function_id) REFERENCES apn_function(id),
    FOREIGN KEY (family_id) REFERENCES apn_family(id)
);

CREATE TABLE IF NOT EXISTS apn_invariant (
    function_id INTEGER PRIMARY KEY,
    walsh_spectrum TEXT,
    delta_rank INTEGER,
    gamma_rank INTEGER,
    multiplier_group_order INTEGER,
    FOREIGN KEY (function_id) REFERENCES apn_function(id)
);




INSERT OR IGNORE INTO apn_family (id, name, description) VALUES
    (1, 'Gold', 'Gold power family'),
    (2, 'Kasami', 'Kasami power family'),
    (3, 'Welch', 'Welch power family'),
    (4, 'Niho', 'Niho power family'),
    (5, 'Inverse', 'Inverse power family'),
    (6, 'Dobbertin', 'Dobbertin power family'),
    (7, 'C1-C2', 'Quadratic APN family C1-C2'),
    (8, 'C3', 'Quadratic APN family C3'),
    (9, 'C4', 'Quadratic APN family C4'),
    (10, 'C5', 'Quadratic APN family C5'),
    (11, 'C6', 'Quadratic APN family C6'),
    (12, 'C7-C9', 'Quadratic APN family C7-C9'),
    (13, 'C10', 'Quadratic APN family C10'),
    (14, 'C11', 'Quadratic APN family C11');
