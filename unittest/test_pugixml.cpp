
#include <gmock/gmock.h>
#include <pugixml.hpp>

using namespace testing;

TEST(TEST_pugixml, DOC) {
    pugi::xml_document doc;

    auto decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version") = "1.0";
    decl.append_attribute("encoding") = "UTF-8";


    auto foo_node = doc.append_child("foo");
    foo_node.append_attribute("bar") = "baz";
    auto call_node = foo_node.append_child("call");

    call_node.append_child(pugi::node_pcdata).set_value("hey");


    std::cout << std::endl;
    doc.save(std::cout);
    std::cout << std::endl;
}






