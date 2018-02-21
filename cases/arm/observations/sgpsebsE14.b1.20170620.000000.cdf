CDF   0   
      time             command_line      sebs_ingest -s sgp -f E14      process_version       ingest-sebs-1.5-0.el6      ingest_software       ingest-sebs-1.5-0.el6      dod_version       sebs-b1-1.4    site_id       sgp    facility_id       E14: Lamont, Oklahoma      
data_level        b1     input_source      ?/data/collection/sgp/sgpsebsE14.00/Table30.20170620_000000.raw     resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

resolution for lat = 0.001
resolution for lon = 0.001
resolution for alt = 1     averaging_interval        30 minutes     sampling_interval         
5 seconds      serial_number         N/A    cdl_program_signature         43690
     qc_standards_version      1.0    	qc_method         Standard Mentor QC     
qc_comment       The QC field values are a bit packed representation of true/false values for the tests that may have been performed. A QC value of zero means that none of the tests performed on the value failed.

The QC field values make use of the internal binary format to store the results of the individual QC tests. This allows the representation of multiple QC states in a single value. If the test associated with a particular bit fails the bit is turned on. Turning on the bit equates to adding the integer value of the failed test to the current value of the field. The QC field's value can be interpretted by applying bit logic using bitwise operators, or by examing the QC value's integer representation. A QC field's integer representation is the sum of the individual integer values of the failed tests. The bit and integer equalivalents for the first 5 bits are listed below:

bit_1 = 00000001 = 0x01 = 2^0 = 1
bit_2 = 00000010 = 0x02 = 2^1 = 2
bit_3 = 00000100 = 0x04 = 2^2 = 4
bit_4 = 00001000 = 0x08 = 2^3 = 8
bit_5 = 00010000 = 0x10 = 2^4 = 16      qc_bit_1_description      !Value is equal to missing_value.       qc_bit_1_assessment       Bad    qc_bit_2_description      "Value is less than the valid_min.      qc_bit_2_assessment       Bad    qc_bit_3_description      %Value is greater than the valid_max.       qc_bit_3_assessment       Bad    qc_bit_4_description      DDifference between current and previous values exceeds valid_delta.    qc_bit_4_assessment       Indeterminate      
datastream        sgpsebsE14.b1      history       Zcreated by user dsmgr on machine ruby at 2017-06-20 01:21:01, using ingest-sebs-1.5-0.el6         G   	base_time                string        2017-06-20 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         Pt   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2017-06-20 00:00:00 0:00          P�   time                	long_name         Time offset from midnight      units         'seconds since 2017-06-20 00:00:00 0:00          P�   qc_time                 	long_name         :Quality check results on field: Time offset from midnight      units         	unitless       description       vThis field contains bit packed values which should be interpreted as listed. No bits set (zero) represents good data.      bit_1_description         9Delta time between current and previous samples is zero.       bit_1_assessment      Indeterminate      bit_2_description         fDelta time between current and previous samples is less than the delta_t_lower_limit field attribute.      bit_2_assessment      Indeterminate      bit_3_description         iDelta time between current and previous samples is greater than the delta_t_upper_limit field attribute.       bit_3_assessment      Indeterminate      delta_t_lower_limit       @�        delta_t_upper_limit       @�,        prior_sample_flag               comment       �If the 'prior_sample_flag' is set the first sample time from a new raw file will be compared against the time just previous to it in the stored data. If it is not set the qc_time value for the first sample will be set to 0.         P�   down_short_hemisp                   	long_name         -Downwelling shortwave hemispheric irradiance       units         W/m^2      	valid_min                	valid_max         D�     
resolution        =���   missing_value         �<         P�   qc_down_short_hemisp                	long_name         MQuality check results on field: Downwelling shortwave hemispheric irradiance       units         	unitless       description       7See global attributes for individual bit descriptions.          P�   up_short_hemisp                 	long_name         +Upwelling shortwave hemispheric irradiance     units         W/m^2      	valid_min                	valid_max         D�     
resolution        =���   missing_value         �<         P�   qc_up_short_hemisp                  	long_name         KQuality check results on field: Upwelling shortwave hemispheric irradiance     units         	unitless       description       7See global attributes for individual bit descriptions.          P�   	down_long                   	long_name         Sky longwave irradiance    units         W/m^2      	valid_min                	valid_max         DH     
resolution        =���   missing_value         �<         P�   qc_down_long                	long_name         8Quality check results on field: Sky longwave irradiance    units         	unitless       description       7See global attributes for individual bit descriptions.          P�   up_long                 	long_name         Surface longwave irradiance    units         w/m^2      	valid_min                	valid_max         DH     
resolution        =���   missing_value         �<         P�   
qc_up_long                  	long_name         <Quality check results on field: Surface longwave irradiance    units         	unitless       description       7See global attributes for individual bit descriptions.          P�   surface_soil_heat_flux_1                	long_name         Surface soil heat flux 1       units         W/m^2      
resolution        =���   missing_value         �<    	valid_min         �H     	valid_max         B�          P�   qc_surface_soil_heat_flux_1                 	long_name         9Quality check results on field: Surface soil heat flux 1       units         	unitless       description       7See global attributes for individual bit descriptions.          P�   surface_soil_heat_flux_2                	long_name         Surface soil heat flux 2       units         W/m^2      
resolution        =���   missing_value         �<    	valid_min         �H     	valid_max         B�          P�   qc_surface_soil_heat_flux_2                 	long_name         9Quality check results on field: Surface soil heat flux 2       units         	unitless       description       7See global attributes for individual bit descriptions.          P�   surface_soil_heat_flux_3                	long_name         Surface soil heat flux 3       units         W/m^2      	valid_min         �H     	valid_max         B�     
resolution        =���   missing_value         �<         P�   qc_surface_soil_heat_flux_3                 	long_name         9Quality check results on field: Surface soil heat flux 3       units         	unitless       description       7See global attributes for individual bit descriptions.          P�   soil_moisture_1                 	long_name         Soil moisture 1, gravimetric       units         %      
resolution        =���   missing_value         �<    	valid_min                	valid_max         BH          P�   qc_soil_moisture_1                  	long_name         =Quality check results on field: Soil moisture 1, gravimetric       units         	unitless       description       7See global attributes for individual bit descriptions.          P�   soil_moisture_2                 	long_name         Soil moisture 2, gravimetric       units         %      
resolution        =���   missing_value         �<    	valid_min                	valid_max         BH          P�   qc_soil_moisture_2                  	long_name         =Quality check results on field: Soil moisture 2, gravimetric       units         	unitless       description       7See global attributes for individual bit descriptions.          P�   soil_moisture_3                 	long_name         Soil moisture 3, gravimetric       units         %      
resolution        =���   missing_value         �<    	valid_min                	valid_max         BH          P�   qc_soil_moisture_3                  	long_name         =Quality check results on field: Soil moisture 3, gravimetric       units         	unitless       description       7See global attributes for individual bit descriptions.          P�   soil_temp_1                 	long_name         Soil temperature 1     units         degC       
resolution        =���   missing_value         �<    	valid_min         �      	valid_max         BH          P�   qc_soil_temp_1                  	long_name         3Quality check results on field: Soil temperature 1     units         	unitless       description       7See global attributes for individual bit descriptions.          P�   soil_temp_2                 	long_name         Soil temperature 2     units         degC       
resolution        =���   missing_value         �<    	valid_min         �      	valid_max         BH          P�   qc_soil_temp_2                  	long_name         3Quality check results on field: Soil temperature 2     units         	unitless       description       7See global attributes for individual bit descriptions.          P�   soil_temp_3                 	long_name         Soil temperature 3     units         degC       
resolution        =���   missing_value         �<    	valid_min         �      	valid_max         BH          P�   qc_soil_temp_3                  	long_name         3Quality check results on field: Soil temperature 3     units         	unitless       description       7See global attributes for individual bit descriptions.          P�   soil_heat_flow_1                	long_name         Soil heat flow 1       units         W/m^2      
resolution        =���   missing_value         �<    	valid_min         �H     	valid_max         B�          Q    qc_soil_heat_flow_1                 	long_name         1Quality check results on field: Soil heat flow 1       units         	unitless       description       7See global attributes for individual bit descriptions.          Q   soil_heat_flow_2                	long_name         Soil heat flow 2       units         W/m^2      
resolution        =���   missing_value         �<    	valid_min         �H     	valid_max         B�          Q   qc_soil_heat_flow_2                 	long_name         1Quality check results on field: Soil heat flow 2       units         	unitless       description       7See global attributes for individual bit descriptions.          Q   soil_heat_flow_3                	long_name         Soil heat flow 3       units         W/m^2      
resolution        =���   missing_value         �<    	valid_min         �H     	valid_max         B�          Q   qc_soil_heat_flow_3                 	long_name         1Quality check results on field: Soil heat flow 3       units         	unitless       description       7See global attributes for individual bit descriptions.          Q   corr_soil_heat_flow_1                   	long_name         .Soil heat flow 1, corrected for soil moisture      units         W/m^2      
resolution        =���   missing_value         �<         Q   qc_corr_soil_heat_flow_1                	long_name         NQuality check results on field: Soil heat flow 1, corrected for soil moisture      units         	unitless       description       7See global attributes for individual bit descriptions.          Q   corr_soil_heat_flow_2                   	long_name         .Soil heat flow 2, corrected for soil moisture      units         W/m^2      
resolution        =���   missing_value         �<         Q    qc_corr_soil_heat_flow_2                	long_name         NQuality check results on field: Soil heat flow 2, corrected for soil moisture      units         	unitless       description       7See global attributes for individual bit descriptions.          Q$   corr_soil_heat_flow_3                   	long_name         .Soil heat flow 3, corrected for soil moisture      units         W/m^2      
resolution        =���   missing_value         �<         Q(   qc_corr_soil_heat_flow_3                	long_name         NQuality check results on field: Soil heat flow 3, corrected for soil moisture      units         	unitless       description       7See global attributes for individual bit descriptions.          Q,   soil_heat_capacity_1                	long_name         Soil heat capacity 1       units         MJ/m^3/degC    
resolution        =���   missing_value         �<         Q0   qc_soil_heat_capacity_1                 	long_name         5Quality check results on field: Soil heat capacity 1       units         	unitless       description       7See global attributes for individual bit descriptions.          Q4   soil_heat_capacity_2                	long_name         Soil heat capacity 2       units         MJ/m^3/degC    
resolution        =���   missing_value         �<         Q8   qc_soil_heat_capacity_2                 	long_name         5Quality check results on field: Soil heat capacity 2       units         	unitless       description       7See global attributes for individual bit descriptions.          Q<   soil_heat_capacity_3                	long_name         Soil heat capacity 3       units         MJ/m^3/degC    
resolution        =���   missing_value         �<         Q@   qc_soil_heat_capacity_3                 	long_name         5Quality check results on field: Soil heat capacity 3       units         	unitless       description       7See global attributes for individual bit descriptions.          QD   energy_storage_change_1                 	long_name         .Change in energy storage 1, 0-5 cm soil layer      units         W/m^2      
resolution        =���   missing_value         �<         QH   qc_energy_storage_change_1                  	long_name         NQuality check results on field: Change in energy storage 1, 0-5 cm soil layer      units         	unitless       description       7See global attributes for individual bit descriptions.          QL   energy_storage_change_2                 	long_name         .Change in energy storage 2, 0-5 cm soil layer      units         W/m^2      
resolution        =���   missing_value         �<         QP   qc_energy_storage_change_2                  	long_name         NQuality check results on field: Change in energy storage 2, 0-5 cm soil layer      units         	unitless       description       7See global attributes for individual bit descriptions.          QT   energy_storage_change_3                 	long_name         .Change in energy storage 3, 0-5 cm soil layer      units         W/m^2      
resolution        =���   missing_value         �<         QX   qc_energy_storage_change_3                  	long_name         NQuality check results on field: Change in energy storage 3, 0-5 cm soil layer      units         	unitless       description       7See global attributes for individual bit descriptions.          Q\   albedo                  	long_name         Albedo     units         	fraction       	valid_min                	valid_max         ?�     
resolution        <#�
   missing_value         �<         Q`   	qc_albedo                   	long_name         'Quality check results on field: Albedo     units         	unitless       description       7See global attributes for individual bit descriptions.          Qd   net_radiation                   	long_name         Net radiation      units         W/m^2      
resolution        =���   missing_value         �<    	valid_min         �H     	valid_max         Dz          Qh   qc_net_radiation                	long_name         .Quality check results on field: Net radiation      units         	unitless       description       7See global attributes for individual bit descriptions.          Ql   surface_soil_heat_flux_avg                  	long_name         .Surface soil heat flux, average of fluxes 1-3      units         W/m^2      
resolution        =���   missing_value         �<    	valid_min         �H     	valid_max         B�          Qp   qc_surface_soil_heat_flux_avg                   	long_name         NQuality check results on field: Surface soil heat flux, average of fluxes 1-3      units         	unitless       description       7See global attributes for individual bit descriptions.          Qt   surface_energy_balance                  	long_name         Surface energy balance     units         W/m^2      	valid_min         ��     	valid_max         D�     
resolution        =���   missing_value         �<         Qx   qc_surface_energy_balance                   	long_name         7Quality check results on field: Surface energy balance     units         	unitless       description       7See global attributes for individual bit descriptions.          Q|   wetness                 	long_name         Wetness, rain detector     units         V      	valid_min         ?Y��   	valid_max         @Fff   
resolution        <#�
   missing_value         �<    comment       B3 V indicates sensor is dry, 1 V indicates sensor is fully wetted           Q�   
qc_wetness                  	long_name         7Quality check results on field: Wetness, rain detector     units         	unitless       description       7See global attributes for individual bit descriptions.          Q�   temp_net_radiometer                 	long_name         Net radiometer temperature     units         degC       	valid_min         �      	valid_max         BH     
resolution        <#�
   missing_value         �<         Q�   qc_temp_net_radiometer                  	long_name         ;Quality check results on field: Net radiometer temperature     units         	unitless       description       7See global attributes for individual bit descriptions.          Q�   battery_voltage                 	long_name         Battery voltage    units         V      	valid_min                	valid_max         Ap     
resolution        <#�
   missing_value         �<         Q�   qc_battery_voltage                  	long_name         0Quality check results on field: Battery voltage    units         	unitless       description       7See global attributes for individual bit descriptions.          Q�   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�          Px   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4          P|   alt              	long_name         Altitude above mean sea level      units         m           P�YHe�Bm�����C��                     C���    B��|    C�Z�    C�w+    ���    �F�    ��H�    A�49    A���    A���    A��{    A��    A��    ���}    �D�    ���?    ���E    �F \    �ÇU    @|[    ?�l    @�    ;t�    :o��    ���    >\
�    C#3    �U�    Cۦ    @;��    B�$    AO+�    @�      @�          Ca��    BL�&    C��    C�    ��V�    ���    ���g    A���    A��<    A��'    A��    A��    A���    ���V    ���    ����    ��э    ��@    ��iY    @:�    ?�m]    @�    ��0�    >C�\    >�>`    >g�`    B��v    ��jj    B�xl    @<4n    A�>w    AP��    @�      @�          C��    B��    Cá    C�y    �J!�    ��v    �U]    A�2-    A���    AƐ�    A��    A��    A�ӏ    �FV�    ���    �N��    �Ib�    ��    �Ry�    @	�    ?�p�    @�    �>��    ��q    �'!\    >x�    A�\)    ��    Aj�    @9�a    A�.    AOB[    @�     @�         BFu    AV_p    C��    C�J    ���    ����    ����    A��X    A��3    A�H�    A�f�    A���    A�o�    �d�    ��/E    �̣    ���    ��P	    �M    @��    ?�s�    @�    ?S�    >�    ?���    >���    ��O    �,�.    �$��    @1��    A�0�    AP�r    @�      @�          ?y��    >>��    C�x�    C�      ����    ��    ����    A�hs    A�ʌ    A�	    A��    A�{�    A��	    ��Q/    �96    �)��    ���B    �:0+    �,�>    @��    ?�y}    @��    ��Bq    >�ة    >��    >D��    ~    ��F�    Y    @"[�    A��M    AO�    @��     @��         ��>   ��w�   C���    Cژ�    �V��    ?j[    @)�m    A�	    A�Υ    A��]    A���    A�I    A�,    ���D    >��    ?�Б    ���u    >�X�    ?��    @�z    ?�}    @�    ?B�[    >��    ?��    >�̅    ��    ?h�0    �    ?�uO    AȬ�    AP�M    @�     @�         �̬G   ���
   C���    C׉�    �uh�    @~�    @|��    A�Ţ    A�ҽ    Aż6    A��~    A�y>    A�3�    �d�    @G�~    @R    �h*    @I*0    @"
(    @gw    ?進    @d�    �Vi    ?U��    ?�Uq    >aܔ    W?    @�/    DM    ?�Q�    A���    AQ�    @Ȝ     @Ȝ         �Ǚ   ��U�   C�°    C��D    ?�t    @��(    @�4�    A�`�    A��    AŹ�    A�,    Aå�    A��	    ��2    @�%F    @J�v    �ԙ    @�1�    @Nr2    @:�    ?�D    @c^    ?Bʛ    ?^�    ?���    >`�	    �|��    @qQ/    �n    ?�7    A���    AO�j    @�      @�          ��2M   ��ޔ   C��    C�0�    ?Y�    @��    @�v`    A�    A��p    AŞ�    A��     A���    A���    ?0�    @Ŕ    @[    ?2��    @���    @_$     @�    ?��    @W�    >�n    ?)Wy    ?��?    >{О    �u��    @��O    �d|�    ?��#    A��*    AP��    @Ϥ     @Ϥ         ����   ��A5   C�YX    C�O\    @��    A    @ĀI    A���    A��    A�~�    A�:�    A�>w    AC    ?�֡    @��'    @��9    ?�v    @�S    @�r�    @��    ?��    @IR    ?G~N    ?�L�    ?�6    >H�Z    �y�    @��N    �cQ     ?��    A��    AP��    @є     @є         ��W   ����   C�R    C��y    @&T�    A�R    @�u�    A��    A��y    A�c�    A���    A��{    A�`    @ )     @��    @�W�    @k    A �7    @��q    @��    ?锯    @=2    ?�    ?!�    ?��p    >u�    �~�^    @�˧    �fQ    ?���    A�ʌ    AO��    @�V     @�V         ���M   ����   C�*�    C���    @K�    AqA    @���    A�>w    A��    A�L0    A�N�    A��z    A�Dg    @$��    A��    @�m�    @'�    A�;    @��+    @�x    ?�=    @2�    ?4�    ?^@d    ?~�A    >Tm�    �w>�    @�hI    �\��    ?��W    A�)�    AOxl    @�     @�         ���h   ��1Q   C���    CѸ�    @��    A&-�    Agm    A�\    A��l    A��    A���    A��T    A�6�    @S6    A�    @�G�    @V	    A\    @�,|    @�    ?�    @v    ?F̒    ?1�    ?։"    >W�$    �w�	    @�    �YtT    ?��    A��    APp�    @��     @��         ���   ��d   C��=    CИ1    @�D|    A5��    Ał    A��t    A�      A���    A��    A�&�    A�@    @��    A'�    @�in    @��r    A(�    @�F    @|�    ?�     @�    ?B0�    ?^Y�    ?�ڥ    >Wƹ    �~�)    A�s    �]p    ?��    A��V    APp�    @؜     @؜         ��!�   ��    C��F    C��=    @�2v    A<A�    A	զ    A�l�    A�    A���    A�>w    A�e�    A��    @� �    A0ƨ    @�a    @�DR    A1�K    @���    @\>    ?鬛    @
��    ?���    ?%u)    ?��F    >_��    �y*K    A�    �V'    ?��    A��?    AN�A    @�^     @�^         ���   �|Z>   C�C    C�/�    @�>�    AJ    A�U    A� �    A�1    A���    A���    A���    A�у    @�m3    A7��    @�Bp    @���    A9 �    @�    @:~    ?�6    @
�    ?<�    ?��    ?�`W    >G�    �w�)    A0    �S^�    ?�    A�U�    AO�E    @�      @�          ��j   ���#   C�aH    C���    @���    AS�R    AE�    A��]    A�    Aģ�    A��    A���    A�ӏ    @���    AD��    @�"�    @�LD    AE�}    @�Ɇ    @    ?�}    @
��    ?��l    ?^{�    ?��    >^n�    �}Q�    A~    �WR�    ?��    A��"    AOĜ    @��     @��         ��Ë   �5��   C�I�    C̕`    @�f'    A\8    A0U    A��$    A��    AĦL    A�'R    A���    A��N    @��    APA     A�O    @�{t    AQ��    A��    @�7    ?鷿    @
��    >X     ?)�P    ?���    >��    �{^�    A�"    �Sf�    ?��%    A��    APy�    @ߤ     @ߤ         ���   �MZ\   C�T�    C�q    A
    Aj+    A9�    A�O    A��    A�nc    A�r�    A�5?    A��    @�҉    AW�a    A�     @Ͽ�    AY(�    Ah�    @�D    ?�^    @
�    ?��6    ?�}�    ?��'    >+Zd    �xΥ    A-$�    �NK�    ?�8�    A���    AP�K    @�     @�         ��{�   ����   C���    C�.�    @��    As_    A+�    A�J    A�_    A�^5    A��t    A�O�    A�֡    @�Mj    Aa=    A
v�    @�z�    Ab�-    A^    @��    ?龡    @
��    ?�)    ?�f�    ?��[    >�3�    �t�    A2҉    �GQ�    ?�)�    A��    AOj�    @�     @�         ���S   �/�l   C���    Cˎ�    Av�    Au��    A��    A�Ɇ    A�	    A�^�    A���    A�v`    A��y    @꿜    Af[W    A	��    @��    Ag�    A%    @��    ?���    @
�2    ?v�,    ?^i8    ??h    >Ë    �o�    A0�e    �B��    ?�i�    A���    AP5�    @�u     @�u         ����   ��d�   C�iX    C��    A��    Ar�w    A�    A��'    A�    A�=    A�x    A��'    A��    @���    Af_    A	��    @��    Ag�/    A9    @�q    ?���    @
�    ?z��    ?.�    ?���    >[nP    �s�f    A3"�    �G    ?�O�    A�V�    AP�    @�V     @�V         ��s   ��(�   C���    C�}q    @���    Ar��    A~]    A�r�    A� �    A�$    A�Gz    A��2    A�3�    @���    Af]d    As    @�M    Agی    A	��    @{     ?��*    @
�%    ?u!    ?2=�    ?�d�    >W��    �s��    A/)�    �Gj�    ?��    A�u%    APy�    @�7     @�7         @��0    ?��    C�@b    C�+#    A	u�    Ax��    A�    A�D3    A�#    A�-    A���    A�9�    A�[�    @��    AiE9    A	8R    @�V    Ajȴ    A��    @f{    ?��!    @
�    ?�r�    ?^,    ?T�    >E�    �I�-    A3��    �r-    ?�=    A���    AP8    @�     @�         B�O�    A9c�    C�K�    C�D�    @�m�    Ar�B    A�U    A�$@    A�&    A�'�    A�՛    A�yr    A���    @��    Ah�    A�I    @��    Aj7�    AV    @Xd    ?���    @
��    >�f    ?	y}    ?X4	    >8�:    ��*    A,J�    @�=�    ?��    A��    AO($    @��     @��         C��    A�=    C�2�    Cο    @��    Ad�\    A	��    A�    A�'R    A�6z    A�o     A��    A�$�    @��_    AWJ�    @��D    @��    AX�!    @�p;    @H�    ?�˼    @
�,    ?Ϣ    ?=�    ?h��    >Dr;    Blյ    A#    B��=    ?�u    A��.    AO�y    @��     @��         C|~�    BGy�    C��    C�>�    @��    ADJ#    @���    A���    A�,    A�o    A���    A�p    A���    @ǎ�    A9�H    @��+    @�Y6    A:�    @֢^    @G0    ?���    @
�H    ?U�3    ?_5    =.��    >J@S    C�    A�	    C�Z    ?���    A�K^    APc    @�     @�         C�@�    B��    C�{d    C۩7    @��,    A�?    @�$_    A��2    A�.�    A��    A�m�    A�&L    A��q    @�o�    AP3    @��H    @��M    AJX    @�1Q    @5+    ?�Ҟ    @
�<    >�s    >a�@    =y�@    >F�    CR��    @⁮    CY�b    ?�*�    A��	    AP�2    @�     @�         C�d    B�n    C�_�    C��b    @�Vm    @�Ln    @{�^    A��0    A�5    A��K    A�1�    A�G    A��    @�)5    @���    @��    @�X    @�O�    @�k�    @!�    ?��0    @
�l    >gD,    <��
    �.7X    >B½    C�u    @�4Y    C��    ?���    A��[    AP?    @�}     @�}         D�L    B�=�    C��    C��    @X�    @���    @N�    A�]d    A�49    A��[    A��    A��    A��i    @R�A    @�5    @#G�    @U�	    @�2�    @&F�    @ ��    ?��4    @
��    =L�A    ��Z�    ����    >=��    C�1H    @W �    C��    ?��    A�&�    AP��    @�^     @�^         D"�5    B�
    C�֨    C�l)    ?�4    @q    <֔a    A�!�    A�0�    AÔ    A�O    A��    A�=�    @Ց    @?�    ?P��    @��    @�    ?T��    @ �`    ?��A    @
o    ��.U    = ��    �N%    >9;�    C�s    ?��    Cܬ�    @pe    A��c    AP�    @�?     @�?         D7d    C?�    C��    C�\J    >��k    �Ȯ�    ��    A���    A�,�    A�x�    A�v+    A�Ta    A��    ?0��    ��Q�    ��Z    ?3ew    ��2v    ��]y    @ ��    ?�Ц    @
b�    ���    � �    ��n�    >5�    C��    ��D(    C�v�    @1e    A�.    AO{�    @�      @�          DJ@    C�F    C�=    C�F    ���    ��y�    ���6    A��y    A�(�    A�h�    A�ʌ    A��s    A�u    �PC    ��l"    �iT�    �R��    ��b�    �m�J    @ ��    ?��    @
[�    �!>    ��k�    �Ϗ�    >2n�    D��    �y��    D�B    @7U�    A���    AQC    @�     @�         DZB    C'�    C��o    C� �    �:��    ��j    ����    A��t    A�%�    A�Gz    A�`�    A�b    A�W�    �o    �m�    ���d    �
M+    �[�    ���h    @ �5    ?��m    @
L�    �A�i    �Y�X    ����    >1M�    DZ�    ��H,    D�V    @:3    B �C    AO��    @��     @��         Dfu    C��    C��J    D    ���r    �}b�    � q    A��x    A� �    A� �    A�1    A�?H    A��<    �x�    �k?    ����    �{}�    �l�m    ��2�    @ �    ?��*    @
;�    �9�    ���    �^t    >.��    D!�    �"�U    D�    @;rq    Bl"    AP��    @��     @��         Dl�X    C��    C���    D��    ��d�    ��T�    �Al�    A�q    A��    A���    A��    A�.�    A��b    ��H�    ��(    �&    ���a    ��%�    � ��    @ ��    ?���    @
#    ��4    ���    ��o    >,�b    D&��    �]�    D#Sd    @<�    B��    AP�    @�     @�         Dt�    C$    C���    D5    ��E    ��ר    �j;�    A�OB    A��    A�    A�M    A�hs    A�d&    ���    ���    �?�I    ��_[    ��]d    �CM    @ ��    ?�    @	�    ����    �Ǥ@    ���    >+r�    D-�'    ���6    D)�    @;�7    B	�    AP�^    @�B�    @�B�        Dv��    C%d    C�P!    D9'    �$֡    ��4�    ���	    A��    A��    A�O    A�2�    A���    A�X�    ���    �ӿH    �e��    ��N    ��~    �i�    @ r�    ?�    @	�j    ���D    ��rG    ���    >+��    D.1X    ���B    D)�     @<9    B�    AP�,    @�     @�         DmT    C.�    C��`    D�    �@Xy    ����    ����    A�	�    A�	�    A�!b    A�h�    A�.�    A�c�    �%��    ���    ����    �(�    ��"h    ���    @ h�    ?鱅    @	�    ��&�    �ǣ�    ���    >+�T    D&K    ��%�    D!!�    @<��    B1�    AO��    @�#�    @�#�        Dt(�    C$��    C��N    D�1    �T�W    � �    ��}�    A��p    A��    A��    A��    A���    A�H    �6�]    ���    ��v+    �9S�    ��M    ���    @ V    ?�K    @	�s    ���L    �Ǡ�    ��    >,��    D+F�    ��>B    D%�    @<�I    Bn�    AP_p    @�     @�         Dl	'    C�9    C��=    D��    �Z��    �2�    ���    A���    A���    A�p�    A�(�    A�@    A���    �@1'    ����    ���    �B�m    ���f    ��zx    @ 1�    ?�s    @	{t    ��4n    ��֌    ��    >-5�    D$'
    ����    D��    @<�    B��    AP	    @��    @��        D`\�    C�5    CӪ=    D�5    �g��    ��+�    ��Vm    A�f�    A���    A�(X    A���    A���    A�s�    �I��    ��6    ���J    �L��    ��7    ��u�    @  G    ?�9    @	[W    ��)     ��J#    ��_    >.�m    D2�    ��.    D�\    @<�    BX�    AP^�    @�u     @�u         DS>�    C�    C���    D    �j˒    ����    ���
    A�-    A��"    A��    A���    A�P    A��    �O��    �� 4    ����    �R�\    ���    ��6�    @ �    ?�N    @	;�    ���    ��^J    ���    >1�t    D    ��M6    D
�    @<�l    B�    AP]�    @��    @��        DB��    C	��    Cӑ�    Dnf    �_n�    ����    ��A�    A���    A��8    A��{    A�*e    A�^    A�P}    �K+�    ��Em    �zO    �M�j    �۫�    �~�    ?�خ    ?钸    @	    ���    ��^J    ��(�    >5�    D�    ����    C��    @<�    BX�    AP�K    @�V     @�V         D0)�    C VF    C�r-    D�P    �W��    �ϡ-    �qTa    A���    A��    A�=<    A�<�    A�n�    Aȇ�    �@:�    ��~    �V˒    �B�    ��Q    �Z�S    ?��<    ?�u    @��    ��5    ���    ���    >:}$    C��f    ��k    Cܷ
    @<��    B	    AP��    @�ƀ    @�ƀ        DӦ    B�.�    C��}    D;    �;�4    ��lW    �"e�    A�q�    A��/    A��-    A�$    A�<j    A�O�    �)�W    ��5�    �,=    �,(�    ��Dg    �ܼ    ?�f�    ?��    @��    �uu\    �     ��+=    >?�-    C��    �f?�    C��    @;��    B��    APƨ    @�7     @�7         D��    B�W
    Cё�    D��    �P�    ���    �ֲ�    A�Y    A��    A�p�    A���    AŨ$    Aɀi    �	n    ��S�    ��
g    �:    ��%�    �Р�    ?�]    ?醘    @��    �Al�    ��C|    �BD�    >G�    C��    �.*�    C�&�    @;b9    BJ=    AP�1    @���    @���        C�+d    B��+    CϮ�    D �)    ��<`    �E�    �c��    A���    A��s    A��    A�P    A�ʌ    A�YK    ��Ri    �C�    �w��    ��4�    �EO    �{�    ?��    ?��    @i/    �@��    9�    >�(V    >QV    Cr��    ���    Ck;�    @;�T    B�m    AP��    