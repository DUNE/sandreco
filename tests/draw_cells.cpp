    // if (plots) {
    //   h_out->cd();
    //   int p[9] = {100, -2000, 2000, 100, -4000, -0, 100, 22500, 26000};

    //   sand_reco::tracker::DigitCollection::fillMap(digits);
    //   sand_reco::tracker::ClusterCollection clusters(&sand_geo, sand_reco::tracker::DigitCollection::getDigits(), sand_reco::tracker::ClusterCollection::ClusteringMethod::kCellAdjacency);
    //   auto digit_vec =  sand_reco::tracker::DigitCollection::getDigits();
    //   std::string tracker_name = digit_vec.begin()->det;

    //   TCanvas* canvas_cluster = new TCanvas("canvas_cluster","canvas_cluster",2000,1000);
    //   canvas_cluster->Divide(2,1);
      
    //   TH2D* h_cluster_yz = new TH2D("h","h", p[6],p[7], p[8], p[3],p[4], p[5]);
    //   TH2D* h_cluster_xz = new TH2D("h","h", p[6],p[7], p[8], p[0],p[1], p[2]);
    //   canvas_cluster->cd(1);
    //   h_cluster_yz->Draw();
    //   canvas_cluster->cd(2);
    //   h_cluster_xz->Draw();
    //   canvas_cluster->Print("clu.pdf(","pdf");

    //   std::map<double, std::vector<TVectorD>> z_to_tracklets;

    //   int color = 2;
    //   for (const auto& container:clusters.getContainers()) {
    //     int gg = 0;
    //     for (const auto& cluster_in_container:container->getClusters()) {
    //       gg++;
    //       if (color > 9) color = 2;
          
    //       EDEPTree tree;
    //       tree.InizializeFromEdep(*ev, sand_geo.getTGeoManager());
          
    //       std::vector<EDEPTrajectory> primaryTrj;
    //       tree.Filter(std::back_insert_iterator<std::vector<EDEPTrajectory>>(primaryTrj), 
    //         [](const EDEPTrajectory& trj) { return trj.GetParentId() == -1;} );

    //       for (auto trj:primaryTrj) {
            
    //         if (trj.GetTrajectoryPoints().find(string_to_component[tracker_name]) == trj.GetTrajectoryPoints().end()) {
    //           continue;
    //         }
            
    //         for (auto& point : trj.GetTrajectoryPoints().at(string_to_component[tracker_name])) {
    //           TEllipse* pt_yz = new TEllipse(point.GetPosition().Z(), point.GetPosition().Y(), 1);
    //           TEllipse* pt_xz = new TEllipse(point.GetPosition().Z(), point.GetPosition().X(), 1);
    //           pt_yz->SetFillStyle(0);
    //           pt_yz->SetLineWidth(1);
    //           pt_yz->SetLineColor(1);
    //           pt_xz->SetFillStyle(0);
    //           pt_xz->SetLineWidth(1);
    //           pt_xz->SetLineColor(1);
    //           canvas_cluster->cd(1);
    //           pt_yz->Draw();
    //           canvas_cluster->cd(2);
    //           pt_xz->Draw();
    //         }
    //       }

    //       TVector3 first_point;
    //       TVector3 last_point;
    //       double min_z = 10e8;
    //       double max_z = -10e8;
    //       for (uint d = 0; d < cluster_in_container.getDigits().size(); d++) {
    //         auto digit = sand_reco::tracker::DigitCollection::getDigit(cluster_in_container.getDigits()[d]);
            
    //         if (digit.z > max_z) {
    //           max_z = digit.z;
    //           last_point = TVector3(digit.x, digit.y, digit.z);
    //         }
    //         if (digit.z < min_z) {
    //           min_z = digit.z;
    //           first_point = TVector3(digit.x, digit.y, digit.z);
    //         }
    //       }

    //       auto true_tracklet = getTrueTrackletOfCluster(first_point, last_point, TVector3(0,0,0),  TVector3 (0,0,0), cluster_in_container.getZ());
    //       double z_start = cluster_in_container.getZ();
    //       // Draw tracklets
    //       TVector2 start_true_tracklet_yz(z_start, true_tracklet[0].Y());
    //       TVector2 start_true_tracklet_xz(z_start, true_tracklet[0].X());
    //       double zy_end = z_start + 5 * cos(atan(true_tracklet[1].Y() / true_tracklet[1].Z()));
    //       double zx_end = z_start + 5 * cos(atan(true_tracklet[1].X() / true_tracklet[1].Z()));
    //       double y_end = true_tracklet[0].Y() + 5 * sin(atan(true_tracklet[1].Y() / true_tracklet[1].Z()));
    //       double x_end = true_tracklet[0].X() + 5 * sin(atan(true_tracklet[1].X() / true_tracklet[1].Z()));
    //       TVector2 end_true_tracklet_yz(zy_end, y_end);
    //       TVector2 end_true_tracklet_xz(zx_end, x_end);
          
    //       TLine* line_yz_true_tracklet = new TLine(start_true_tracklet_yz.X(), start_true_tracklet_yz.Y(), end_true_tracklet_yz.X(), end_true_tracklet_yz.Y());
    //       TLine* line_xz_true_tracklet = new TLine(start_true_tracklet_xz.X(), start_true_tracklet_xz.Y(), end_true_tracklet_xz.X(), end_true_tracklet_xz.Y());
    //       line_yz_true_tracklet->SetLineColor(color + 1);
    //       line_yz_true_tracklet->SetLineWidth(3);
    //       line_xz_true_tracklet->SetLineColor(color + 1);
    //       line_xz_true_tracklet->SetLineWidth(3);
          
    //       canvas_cluster->cd(1);
    //       line_yz_true_tracklet->Draw();
    //       canvas_cluster->cd(2);
    //       line_xz_true_tracklet->Draw();
        


    //       for (auto digit:digit_vec) {
    //         auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));
    //         auto h = cell->second.getSize().h;
    //         auto w = cell->second.getSize().w;
    //         TBox* box_yz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().Y() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().Y() + w/2.);
    //         TBox* box_xz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().X() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().X() + w/2.);
    //         box_yz->SetFillStyle(0);
    //         box_yz->SetLineColor(1);
    //         box_yz->SetLineWidth(1);
    //         box_xz->SetFillStyle(0);
    //         box_xz->SetLineColor(1);
    //         box_xz->SetLineWidth(1);
    //         canvas_cluster->cd(1);
    //         box_yz->Draw();
    //         canvas_cluster->cd(2);
    //         box_xz->Draw();

    //         TMarker* mark_yz = new TMarker(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().Y(), 5);
    //         mark_yz->SetMarkerColor(1);
    //         mark_yz->SetMarkerSize(0.5);
    //         canvas_cluster->cd(1);
    //         mark_yz->Draw();

    //         TMarker* mark_xz = new TMarker(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().X(), 5);
    //         mark_xz->SetMarkerColor(1);
    //         mark_xz->SetMarkerSize(0.5);
    //         canvas_cluster->cd(2);
    //         mark_xz->Draw();


    //         for (auto& kk:digit.hindex) {
    //           const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
    //           TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
    //           TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
    //           l_yz->SetLineColor(1);
    //           l_xz->SetLineColor(1);
    //           canvas_cluster->cd(1);
    //           l_yz->Draw();
    //           canvas_cluster->cd(2);
    //           l_xz->Draw();
    //         }
    //       }

    //       std::vector<sand_reco::tracker::DigitID> digits_cluster = cluster_in_container.getDigits();

    //       for (uint d = 0; d < digits_cluster.size(); d++) {
    //         canvas_cluster->cd();

    //         auto digit = sand_reco::tracker::DigitCollection::getDigit(digits_cluster[d]);
    //         auto cell = sand_geo.getCellInfo(sand_geometry::tracker::CellID(digit.did));

            
    //         // Draw cells of cluster
    //         auto h = cell->second.getSize().h;
    //         auto w = cell->second.getSize().w;

    //         TBox* box_yz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().Y() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().Y() + w/2.);
    //         TBox* box_xz = new TBox(cell->second.getWire().getCenter().Z() - h/2., cell->second.getWire().getCenter().X() - w/2., cell->second.getWire().getCenter().Z() + h/2., cell->second.getWire().getCenter().X() + w/2.);
    //         box_yz->SetFillStyle(0);
    //         box_yz->SetLineColor(color);
    //         box_yz->SetLineWidth(1);
    //         box_xz->SetFillStyle(0);
    //         box_xz->SetLineColor(color);
    //         box_xz->SetLineWidth(1);
    //         canvas_cluster->cd(1);
    //         box_yz->Draw();
    //         canvas_cluster->cd(2);
    //         box_xz->Draw();
            
    //         // Draw true drift time of digits in cluster
    //         TEllipse* el_yz = new TEllipse(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().Y(), 
    //                               sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
    //         TEllipse* el_xz = new TEllipse(cell->second.getWire().getCenter().Z(), cell->second.getWire().getCenter().X(), 
    //                               sand_reco::stt::wire_radius + cell->second.getDriftVelocity() * digit.drift_time);
    //         el_yz->SetFillStyle(0);
    //         el_yz->SetLineWidth(1);
    //         el_yz->SetLineColor(1);
    //         el_xz->SetFillStyle(0);
    //         el_xz->SetLineWidth(1);
    //         el_xz->SetLineColor(1);
    //         canvas_cluster->cd(1);
    //         el_yz->Draw();
    //         canvas_cluster->cd(2);
    //         el_xz->Draw();

    //         // Draw hit segments for the cluster
    //         for (auto& kk:digit.hindex) {
    //           const TG4HitSegment& hseg = ev->SegmentDetectors[digit.det].at(kk);
    //           TLine* l_yz = new TLine(hseg.Start.Z(), hseg.Start.Y(), hseg.Stop.Z(), hseg.Stop.Y());
    //           TLine* l_xz = new TLine(hseg.Start.Z(), hseg.Start.X(), hseg.Stop.Z(), hseg.Stop.X());
    //           l_yz->SetLineColor(color);
    //           l_xz->SetLineColor(color);
    //           canvas_cluster->cd(1);
    //           l_yz->Draw();
    //           canvas_cluster->cd(2);
    //           l_xz->Draw();
    //         }
    //       }

    //       color++; 
    //       canvas_cluster->Write();
    //       canvas_cluster->Print("clu.pdf","pdf");
    //       canvas_cluster->Clear();

    //       canvas_cluster->Divide(2,1);
    //       canvas_cluster->cd(1);
    //       h_cluster_yz->Draw();
    //       canvas_cluster->cd(2);
    //       h_cluster_xz->Draw();

    //     }
    //   }
    //   canvas_cluster->Print("clu.pdf)","pdf");

    // }