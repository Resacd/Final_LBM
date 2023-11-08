// ImageProcessor.h
#ifndef IMAGEPROCESSOR_H
#define IMAGEPROCESSOR_H

#include <opencv2/opencv.hpp>
#include <vector>

class ImageProcessor {
public:
    ImageProcessor(const std::string& imageFile);
    void displayImage();
    void extractDataFromImage(const std::string& imageFile, double minVFR, double maxVFR,
        double minTime, double maxTime, std::vector<double>& vfrData, std::vector<double>& timeData);
    std::string fitVelocityEquation(const std::vector<double>& vfrData, const std::vector<double>& timeData);

private:
    std::string imageFile;
    std::vector<cv::Point2f> points;

    static void onMouse(int event, int x, int y, int flags, void* userdata);
};

#endif
