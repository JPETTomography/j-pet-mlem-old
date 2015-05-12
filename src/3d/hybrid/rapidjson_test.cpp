#include<iostream>
#include<cstdio>

#include "util/test.h"

#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

TEST("rapid_json") {
  FILE* in=fopen("src/3d/hybrid/point_source.js","r");
  if (!in) {
    std::cerr << "cannot open src/3d/hybrid/point_source.js\n";
  }

  rapidjson::Document doc;
  char readBuffer[256];
  rapidjson::FileReadStream input_stream(in, readBuffer, sizeof(readBuffer));
  doc.ParseStream(input_stream);

  rapidjson::Value& obj = doc[0];

  REQUIRE(obj.HasMember("type"));

}
