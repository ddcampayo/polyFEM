    main_data
      << p.x() << "  "            //  0
      << p.y() << "    "          //  1
      << vit->idx.val() << "  "   //  2
      << vit->vol.val() << "  "   //  3
      << vit->alpha.val() << "  " //  4
      << vit->p.val() << "  "     //  5
      << vit->gradp.val() << "  " //  6 & 7
      << vit->U.val() << "  "     //  8 & 9
      << vit->divU.val() << "  "  //  10
      << vit->laplU.val() << "  " //  11 &12
      << vit->Ustar.val() << "  " //  13 & 14
      << vit->Uold.val() << "  "  //  15 & 16
      << vit->pstar.val() << "  "  //  17
      << endl;
