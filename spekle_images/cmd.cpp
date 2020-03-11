#include "image.h"
#include "builder.h"


#include <ctime>
#include <iostream>
using namespace std;

using namespace speckle;

static void input(config & conf) {
	std::cout << "image rows" << std::endl;
	std::cin >> conf.size.rows;

	std::cout << "image cols" << std::endl;
	std::cin >> conf.size.cols;


	std::cout << "spekle number" << std::endl;
	std::cin >> conf.number;

	std::cout << "spekle radius" << std::endl;
	std::cin >> conf.radius;

	std::cout << "output directory" << std::endl;
	std::cin >> conf.filepath;

	if (conf.filepath.length() <= 1)
		conf.filepath = default_config().filepath;

	std::cout << "output directory: " << conf.filepath << std::endl;
}


static void test_image_builder() {
	image img;

	config conf = default_config();
	input(conf);

	image_builder::build_image(conf, img);
}


int main() {
	srand(time(0));
	test_image_builder();


	cout << "input 0 to exit" << endl;
	int n; cin >> n;
}