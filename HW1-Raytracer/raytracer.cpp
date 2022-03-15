#include <iostream>
#include "parser.h"
#include "ppm.h"


int main(int argc, char* argv[])
{

    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    scene.main_render();

}
