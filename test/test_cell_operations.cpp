#include <cassert>
#include <cstring>
#include <cstdio>

#include "../src/globals.h"


void test_cell_id_to_name() {
    printf("test_cell_id_to_name\n");
    
    char name[5];

    cell_id_to_name(0, name);
    assert(! strcmp("A1", name));

    cell_id_to_name(35, name);
    assert(! strcmp("A8", name));

    cell_id_to_name(18, name);
    assert(! strcmp("C4", name));
}


void test_cell_name_to_id() {
    printf("test_cell_name_to_id\n");
    const char * name1 = "A1";
    assert(cell_name_to_id(name1) == 0);

    const char * name2 = "B3";
    assert(cell_name_to_id(name2) == 8);

    const char * name3 = "A8";
    assert(cell_name_to_id(name3) == 35);
}


void test_adjacencies() {
    printf("test_adjacencies\n");

    char name[5];

    assert(N_ADJ[0] == 2);
    cell_id_to_name(ADJ[0][0], name);
    assert(! strcmp(name, "B1"));
    cell_id_to_name(ADJ[0][1], name);
    assert(! strcmp(name, "A2"));

    CellID id = cell_name_to_id("D4");
    assert(N_ADJ[id] == 6);
    cell_id_to_name(ADJ[id][0], name);
    assert(! strcmp(name, "D3"));
    cell_id_to_name(ADJ[id][1], name);
    assert(! strcmp(name, "C4"));
    cell_id_to_name(ADJ[id][2], name);
    assert(! strcmp(name, "E3"));
    cell_id_to_name(ADJ[id][3], name);
    assert(! strcmp(name, "C5"));
    cell_id_to_name(ADJ[id][4], name);
    assert(! strcmp(name, "E4"));
    cell_id_to_name(ADJ[id][5], name);
    assert(! strcmp(name, "D5"));

}


int main() {
    init();
    test_cell_id_to_name();
    test_cell_name_to_id();
    test_adjacencies();
    printf("All tests passed\n");    
    return 0;
}
