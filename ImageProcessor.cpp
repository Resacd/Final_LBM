// ImageProcessor.cpp
#include "ImageProcessor.h"
#include <iostream>
#include <fstream>

ImageProcessor::ImageProcessor(const std::string& imageFile) : imageFile(imageFile) {}

void ImageProcessor::displayImage() {
    cv::Mat a = cv::imread(imageFile);

    if (a.empty()) {
        std::cerr << "Error: Could not open image." << std::endl;
        return;
    }

    cv::imshow("Image", a);
    cv::setMouseCallback("Image", onMouse, this);
    cv::waitKey(0);

    if (points.size() < 4) {
        std::cerr << "Error: Not enough points selected." << std::endl;
        return;
    }

    float x1 = points[0].x;
    float y1 = points[0].y;
    float x2 = points[1].x;
    float y2 = points[1].y;
    float x3 = points[2].x;
    float y3 = points[2].y;

    float Lx, Ly, Lx0, Ly0;

    std::cout << "Length in x-direction? ";
    std::cin >> Lx;
    std::cout << "Length in y-direction? ";
    std::cin >> Ly;
    std::cout << "Real value of x at origin? ";
    std::cin >> Lx0;
    std::cout << "Real value of y at origin? ";
    std::cin >> Ly0;

    std::vector<cv::Point2f> inputPoints(points.begin() + 3, points.end());

    std::ofstream X3File("X3.txt");
    for (const auto& point : inputPoints) {
        float X3 = (point.x - x1) * Lx / (x2 - x1) + Lx0;
        X3File << X3 << std::endl;
    }
    X3File.close();

    std::ofstream Y3File("Y3.txt");
    for (const auto& point : inputPoints) {
        float Y3 = (point.y - y1) * Ly / (y3 - y1) + Ly0;
        Y3File << Y3 << std::endl;
    }
    Y3File.close();
}

void ImageProcessor::onMouse(int event, int x, int y, int flags, void* userdata) {
    ImageProcessor* processor = static_cast<ImageProcessor*>(userdata);

    if (event == cv::EVENT_LBUTTONDOWN) {
        cv::Point2f point(static_cast<float>(x), static_cast<float>(y));
        processor->points.push_back(point);

        if (processor->points.size() >= 500) {
            cv::destroyWindow("Image");
        }
    }
}


void ImageProcessor::extractDataFromImage(const std::string& imageFile, double minVFR,
    double maxVFR, double minTime, double maxTime, std::vector<double>& vfrData, std::vector<double>& timeData) {
    cv::Mat image = cv::imread(imageFile, cv::IMREAD_GRAYSCALE);

    if (image.empty()) {
        std::cerr << "Error: Could not open image." << std::endl;
        return;
    }

    int width = image.cols;
    int height = image.rows;

    vfrData.clear();
    timeData.clear();

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            // Extract pixel intensity (VFR) and calculate time based on x-coordinate
            double vfr = minVFR + (maxVFR - minVFR) * (static_cast<double>(image.at<uchar>(y, x)) / 255.0);
            double time = minTime + (maxTime - minTime) * (static_cast<double>(x) / (width - 1));

            vfrData.push_back(vfr);
            timeData.push_back(time);
        }
    }
}

std::string ImageProcessor::fitVelocityEquation(const std::vector<double>& vfrData, const std::vector<double>& timeData) {
    int numDataPoints = vfrData.size();

    if (numDataPoints == 0) {
        return "No data to fit.";
    }

    double sumTime = 0.0;
    double sumVFR = 0.0;
    double sumTimeSquared = 0.0;
    double sumTimeVFR = 0.0;

    for (int i = 0; i < numDataPoints; ++i) {
        double time = timeData[i];
        double vfr = vfrData[i];

        sumTime += time;
        sumVFR += vfr;
        sumTimeSquared += time * time;
        sumTimeVFR += time * vfr;
    }

    double slope = (numDataPoints * sumTimeVFR - sumTime * sumVFR) / (numDataPoints * sumTimeSquared - sumTime * sumTime);
    double intercept = (sumVFR - slope * sumTime) / numDataPoints;

    std::string equation = "Velocity = " + std::to_string(intercept) + " + " + std::to_string(slope) + " * Time";

    return equation;
}
