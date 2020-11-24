#include "image.h"
#include "builder.h"
#include "stereo.h"
#include "transform.h"


int eigen_examplex();
int run_transform();



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


static void init() {
	srand(static_cast<unsigned int>(time(0)));
}


static void entrance() {
	int n;

select:

	//ÌבÊ¾
	{
		cout << "1: test_image_builder" << endl;
		cout << "2: test_intrinsic_parameters" << endl;
		cout << "3: eigen_example " << endl;
		cout << "4: test_transform " << endl;
		cout << "5: run_transform " << endl;
		cout << "0: quit " << endl;
	}

	cin >> n;

	switch (n)
	{
	case 1:
		test_image_builder();
		break;
	case 2:
		speckle::test_intrinsic_parameters();
		break;
	case 3:
		eigen_examplex();
		break;

	case 4:
		speckle::test_transform();
		break;

	case 5:
		run_transform();
		break;

	case 0:
		return;

	default:
		break;
	}

	goto select;
}





int main() {
	init();
	entrance();



	cout << "input 0 to exit" << endl;
	int n; cin >> n;
}