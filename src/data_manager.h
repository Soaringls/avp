#ifndef AVP_DATAMANAGER_H_
#define AVP_DATAMANAGER_H_

#include <mutex>
#include <utility>
#include <vector>

#include "data_types.h"
#include "frame.h"
#include "osliding_buffer.h"

namespace qcraft {
namespace localization {
namespace visual {
class DataManager {
 public:
  static DataManager* GetInstance() {
    static DataManager instance;
    return &instance;
  }

  ~DataManager() {
    Reset();
    delete perception_buf_;
    delete lidar_buf_;
    delete gnss_buf_;
    delete position_buf_;
#ifdef VISUALIZATION
    delete raw_img_buf_;
#endif
  }

  DataManager(const DataManager&) = delete;
  DataManager& operator=(const DataManager&) = delete;

  void Reset() {
    gnss_buf_->Clear();
    perception_buf_->Clear();
    position_buf_->Clear();
    lidar_buf_->Clear();
#ifdef VISUALIZATION
    raw_img_buf_->Clear();
#endif
  }

  void Put(const int64_t timestamp, const GnssData& data) {
    std::lock_guard<std::mutex> guard(data_lock_gnss_);
    gnss_buf_->Put(timestamp, data);
  }

#ifdef VISUALIZATION
  void Put(const int64_t timestamp, const ImgData& data) {
    std::lock_guard<std::mutex> guard(data_lock_raw_img_);
    raw_img_buf_->Put(timestamp, data);
  }
#endif

  void Put(const int64_t timestamp, const PerceptionData& data) {
    std::lock_guard<std::mutex> guard(data_lock_perception_);
    perception_buf_->Put(timestamp, data);
  }

  void Put(const int64_t timestamp, const LidarPerceptionData& data) {
    std::lock_guard<std::mutex> guard(data_lock_lidar_);
    lidar_buf_->Put(timestamp, data);
  }

  void Put(const int64_t timestamp, const PositionData& data) {
    std::lock_guard<std::mutex> guard(data_lock_position_);
    position_buf_->Put(timestamp, data);
  }

  Frame* GetLastestFrame() {
    Frame* frame = nullptr;
    if (last_time_ == -1) {
      last_time_ = perception_buf_->StartTime();
    }

    int64_t p_time = perception_buf_->EndTime();

    PerceptionData perception_data;
    {
      std::lock_guard<std::mutex> guard(data_lock_perception_);
      if (!perception_buf_->GetNearestElements(p_time, &perception_data))
        return frame;
    }

    // LidarPerceptionData lidar_data;
    // {
    //   std::lock_guard<std::mutex> guard(data_lock_lidar_);
    //   if (!lidar_buf_->GetNearestElements(p_time, &lidar_data)) return frame;
    // }

#ifdef VISUALIZATION
    ImgData img_data;
    {
      std::lock_guard<std::mutex> guard(data_lock_raw_img_);
      if (!raw_img_buf_->GetNearestElements(p_time, &img_data)) return frame;
    }
#endif

    GnssData gnss_data;
    {
      std::lock_guard<std::mutex> guard(data_lock_gnss_);
      if (!gnss_buf_->GetNearestElements(p_time, &gnss_data)) return frame;
    }

    std::vector<std::pair<int64_t, PositionData>> position_datas;
    {
      std::lock_guard<std::mutex> guard(data_lock_position_);
      position_buf_->GetElements(last_time_, p_time, &position_datas);
      if (position_datas.empty()) return frame;
    }

    frame = new Frame();
    frame->gnss = gnss_data;
    frame->perception = perception_data;
    // frame->lidar_data = lidar_data;
    for (const auto& position : position_datas) {
      frame->positions.emplace_back(position.second);
    }
#ifdef VISUALIZATION
    frame->raw_img = img_data;
#endif
    frame->timestamp = p_time;
    frame->loc.timestamp = p_time;

    last_time_ = p_time;
    return frame;
  }

 private:
  DataManager(/* args */) {
    position_buf_ = new SlidingBuffer<int64_t, PositionData>(100, 15);
    perception_buf_ = new SlidingBuffer<int64_t, PerceptionData>(20, 120);
    lidar_buf_ = new SlidingBuffer<int64_t, LidarPerceptionData>(20, 200);
#ifdef VISUALIZATION
    raw_img_buf_ = new SlidingBuffer<int64_t, ImgData>(20, 120);
#endif
    gnss_buf_ = new SlidingBuffer<int64_t, GnssData>(20, 60);
  }

 private:
  SlidingBuffer<int64_t, PerceptionData>* perception_buf_;
  std::mutex data_lock_perception_;
  SlidingBuffer<int64_t, GnssData>* gnss_buf_;
  std::mutex data_lock_gnss_;
  SlidingBuffer<int64_t, PositionData>* position_buf_;
  std::mutex data_lock_position_;
  SlidingBuffer<int64_t, LidarPerceptionData>* lidar_buf_;
  std::mutex data_lock_lidar_;
#ifdef VISUALIZATION
  SlidingBuffer<int64_t, ImgData>* raw_img_buf_;
  std::mutex data_lock_raw_img_;
#endif
  int64_t last_time_{-1};
};
}  // namespace visual
}  // namespace localization
}  // namespace qcraft

#endif
