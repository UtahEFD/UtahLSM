CDF   0   
      time             command_line      sebs_ingest -s sgp -f E14      process_version       ingest-sebs-1.5-0.el6      ingest_software       ingest-sebs-1.5-0.el6      dod_version       sebs-b1-1.4    site_id       sgp    facility_id       E14: Lamont, Oklahoma      
data_level        b1     input_source      ?/data/collection/sgp/sgpsebsE14.00/Table30.20170601_000000.raw     resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

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
datastream        sgpsebsE14.b1      history       Zcreated by user dsmgr on machine ruby at 2017-06-01 01:21:00, using ingest-sebs-1.5-0.el6         G   	base_time                string        2017-06-01 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         Pt   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2017-06-01 00:00:00 0:00          P�   time                	long_name         Time offset from midnight      units         'seconds since 2017-06-01 00:00:00 0:00          P�   qc_time                 	long_name         :Quality check results on field: Time offset from midnight      units         	unitless       description       vThis field contains bit packed values which should be interpreted as listed. No bits set (zero) represents good data.      bit_1_description         9Delta time between current and previous samples is zero.       bit_1_assessment      Indeterminate      bit_2_description         fDelta time between current and previous samples is less than the delta_t_lower_limit field attribute.      bit_2_assessment      Indeterminate      bit_3_description         iDelta time between current and previous samples is greater than the delta_t_upper_limit field attribute.       bit_3_assessment      Indeterminate      delta_t_lower_limit       @�        delta_t_upper_limit       @�,        prior_sample_flag               comment       �If the 'prior_sample_flag' is set the first sample time from a new raw file will be compared against the time just previous to it in the stored data. If it is not set the qc_time value for the first sample will be set to 0.         P�   down_short_hemisp                   	long_name         -Downwelling shortwave hemispheric irradiance       units         W/m^2      	valid_min                	valid_max         D�     
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
   missing_value         �<         Q�   qc_battery_voltage                  	long_name         0Quality check results on field: Battery voltage    units         	unitless       description       7See global attributes for individual bit descriptions.          Q�   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�          Px   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4          P|   alt              	long_name         Altitude above mean sea level      units         m           P�Y/Y Bm�����C��                     C7
�    B �    C��1    C�}    �.��    �%O�    ��8\    A���    A�H    A���    A���    A��7    A�J�    �-��    �$*0    ���l    �/�    �%    ��֡    ?���    ?�>    ?�&�    ;�η    �p��    ?|�(    >^E=    B�3�    �yS    B�D�    @=�    A�T�    AQw2    @�      @�          C-�#    B&�P    CЮ    C�iy    �9�    ��[-    ���    A�o�    A�E    A��A    A��|    A���    A��    ���    ���w    ��d0    ��    ��Y�    ��h�    ?��    ?�;    ?�a    >4��    >��    ?K��    >u�o    B���    ���B    B��(    @=�    A�HK    AO�*    @�      @�          B�E,    A��0    Cα�    C��}    ���`    �ȭ    �/]y    A�_p    A�E    A�xl    A���    A�9�    A���    ��*�    ��$�    �m�M    ��#    ��>l    �m�3    ?�~    ?�;    ?��    >%W    �7Gz    ?x��    >q�    A���    ���    AQ.I    @<��    A�Ţ    AQj    @�     @�         A�ݘ    @��2    Cϴ    C��    �ï�    ���]    �	0    A�Np    A�F�    A�n�    A���    A���    A��    ���    ��%1    �;.    ��x�    ����    �;2    ?�}k    ?�<�    ?��    >��    ?/��    ?H     >HWm    ��W?    ��n    ��5?    @<f    A�f�    AP�<    @�      @�          ��Z\   �y�   C�o�    Cَ�    ��.�    ��    ���e    A�E�    A�H�    A�l"    A�(�    A�~�    A���    ���0    �4,�    ���    ��$    �5,�    ���    ?�u�    ?�>W    ?��    >�W�    >��    ?j:�    ?j]�    �ԭC    �5�S    ��U�    @
u�    A�&�    AP>B    @��     @��         �=ag   ���W   CǺ=    C��1    �7V�    ����    ���.    A�4�    A�J�    A�i    A��q    A���    A��Q    �T�4    �d��    ����    �V*�    �e�'    ���    ?�f{    ?�@O    ?��    >��+    ?PBc    ?jr�    >�g    ��v�    ����    ��$    @��    A�&L    AP�    @�     @�         �=Q�   ���E   CƧ�    C��    ����    ?���    ?�\    A�$�    A�N    A�c     A��    A�'R    A���    ���T    ?vt    >+0�    ���3    ?w�,    >+3�    ?�Xd    ?�CB    ?���    ?��    ?��.    ?�E�    >�<    ���K    ?��8    �ڭ    @3	    A��b    AP��    @Ȝ     @Ȝ         �d�   �l�   C�B    C�c�    >�ץ    @J��    @5�    A�7    A�P    A�^�    A�($    A�Ta    A�3�    =Ტ    @�    ?�#�    =�o;    @�&    ?�&B    ?�N�    ?�D�    ?���    >���    ?,2�    ?���    >���    ���    @	7�    ��G    @�    A���    AP    @�      @�          �tZq   �_p&   C�d�    C�c�    @��    @�-�    @PӮ    A�V    A�Q    A�V    A�u�    A��q    A�    ?vr�    @s��    @	�    ?xX�    @u*�    @	<    ?�Dg    ?�E�    ?��    ?�)    ?$��    ?�~�    >i��    �+��    @PH�    ���    @v�    A�z�    AO_    @Ϥ     @Ϥ         �z��   �u�   C��     Cј�    @�    @��t    @jEN    A��    A�R     A�L�    A��e    A��    A���    ?��    @��    @-�>    ?��M    @���    @-̎    ?�?}    ?�F�    ?��    ?]_g    ?/݇    ?q�    >{�    �*��    @t�    �ں    @D�    A�v`    AP�H    @є     @є         �{�H   ��]   C�}�    C���    @40U    @ׂ�    @���    A��    A�V    A�Em    A��9    A��    A��    ?߃�    @��    @@�.    ?�<6    @��{    @@��    ?�=�    ?�Jw    ?���    ?�$t    ?��    ?�$J    >��     �'�    @��    ���    ?��    A��6    AP,�    @�V     @�V         �u�   �}92   C��    C�}    @*�u    @�8    @tvK    A�M    A�Y�    A�A�    A��K    A��    A���    ?��7    @�J    @Kԕ    ?�t�    @�O    @K�=    ?�;�    ?�Mj    ?��    ?4Lj    ? O    ?"��    >��1    ��C    @�;d    ��F    ?�?�    A�F�    AO5?    @�     @�         �v�   ��j   C��)    Cϲ-    @Z�    @͖S    @���    A��    A�Z�    A�?}    A�O    A�\)    A��    @    @�j    @Kԕ    @&�    @�w    @K�    ?�:�    ?�N�    ?�޾    ?���    ?\�    ?x��    >�'�    ��C    @�W�    ��    ?���    A�{�    AP͟    @��     @��         �m"�   ���   C��q    C�U?    @O�>    @�    @�5~    A�    A�^�    A�AU    A�\�    A��V    A�چ    @H�    @���    @LA_    @p�    @��    @LA�    ?�;:    ?�R     ?��a    ?`�    ?��    ?x�L    >�F�    ��    @�,    ����    ?��    A�k�    AOj    @؜     @؜         �|�   ���   C��V    C΢o    @D\    @�P    @���    A�+    A�aH    A�>B    A���    A���    A�Y    @"B    @�If    @W��    @#��    @�Q    @W�P    ?�>-    ?�TL    ?���    ?i�    ?��    ?O0�    >��q    �
�    @�    ��n    ?�@    A�	    AP��    @�^     @�^         �v��   ��<�   C�3�    C��9    @\�    @�J�    @��G    A��    A�c     A�>�    A���    A�0�    A�33    @.3�    @���    @_�7    @/��    @��j    @_��    ?�9C    ?�U�    ?��    ?4�    ?,r�    ?x�    >�B|    ��    @�,�    ���>    ?�{    A��    AP�<    @�      @�          ���   ����   C���    C�_}    @j�    @�t    @�Ln    A��    A�ff    A�>B    A�M6    A�~�    A�Y    @<#d    @���    @d��    @=��    @�    @d��    ?�9C    ?�X�    ?���    ?4�m    ?\"�    ?O|$    >�:�    �mC    @��    ����    ?�hs    A��X    AP�<    @��     @��         ��N   ��	�   C�!    C�dZ    @h�    @��0    @�΅    A���    A�iy    A�;�    A��0    A�ܒ    A���    @E�[    @ʱ�    @qG�    @G6P    @��k    @qG�    ?�7L    ?�[�    ?��w    ?�;    ?0    >�P    >�P    �ٚ    @�;�    ����    ?���    A���    AP�.    @ߤ     @ߤ         ��b�   ���E   C���    C�)�    @z    @�T�    @���    A�{    A�m�    A�;0    A��    A� �    A��    @Kd0    @� �    @t��    @L��    @�L�    @t��    ?�:�    ?�_    ?��#    ?497    ?0@}    ?"�'    >�VK    �^    @��     ��    ?�}    A���    AP!    @�     @�         ����   ���\   C���    Cɦf    @�p�    @���    @�'�    A��    A�n/    A�<�    A�I�    A���    A�:�    @Y|    @ؚV    @��Y    @[(�    @��$    @��n    ?�=�    ?�_�    ?��r    ?�qv    ? �z    ?|#S    >��    �	�t    @�*    ���v    ?ǋ�    A��    AP	    @�     @�         �xt~   ���Q   C�nV    C��q    @|�    @�b    @��    A��    A�po    A�;0    A��-    A��    A��u    @g�
    @�    @���    @i��    @��    @���    ?�;�    ?�a�    ?��#    >� �    ? �r    ?Lu�    >���    �s�    @��s    ���]    ?���    A��    AP9X    @�u     @�u         �k��   ��5�   C�i�    CȽ/    @�Xy    A ��    @���    A��    A�rG    A�7�    A��    A���    A��d    @�2�    @�զ    @���    @�.�    @�%�    @���    ?�:?    ?�c^    ?���    ?aL4    ?0�    >�    >�^N    ��    @�o�    ��q    ?�H�    A��#    AO�    @�V     @�V         �v�   ����   C�5    CǤ�    @��q    A{�    @��    A��.    A�v+    A�3�    A�lW    A���    A�2�    @���    @    @��O    @���    @���    @���    ?�6�    ?�f�    ?��A    ?a8    ?06    >�1    >�4�    ��ZQ    @���    ���     ?���    A��w    APk�    @�7     @�7         ?�Vm    �E�   C���    C��    @�2�    A�    @��o    A�o    A�xl    A�5�    A��'    A�5?    A��    @��"    @��    @�2a    @��(    @��    @�2    ?�8�    ?�h�    ?���    ?Ո    ?0	-    ?&�    ��?   ��t�    @���    ����    ?���    A�?}    AO��    @�     @�         A� '    @���    C�3    CƄ9    @���    @��    @�p�    A�    A�}V    A�6�    A�1�    A��G    A��5    @�0    @�g�    @���    @�V    @�S    @���    ?�8�    ?�m3    ?���    ?3�    ?0�    ?%��    >:��    ��3�    @��m    ����    ?��    A�]�    AO��    @��     @��         B�]"    A�    C��    Cɽq    @�|�    @߻�    @�[�    A�    A�~�    A�4�    A���    A�#    A�zD    @{��    @�{     @�d0    @}�T    @߻�    @�c�    ?�<6    ?�n/    ?�Ց    >4�z    7s2�    >�~_    >uds    A�Ö    @��
    BW$    ?�(x    A�E�    AP��    @��     @��         C\<j    BSe�    C�:�    C� �    @���    @�h    @k%F    A�_    A��4    A�5�    A�<�    A��    A���    @g��    @��    @t��    @iVC    @�-�    @t�Y    ?�>-    ?�o�    ?��9    ??L+    ?�    ��    >u��    B�r�    @�d0    C�5    ?��    A��    AO�m    @�     @�         C�`     B���    C�_;    C��)    @R|�    @yN�    @5�,    A�     A���    A�7    A��*    A�l�    A���    @5�g    @h�b    @X    @7Y`    @i��    @    ?�F�    ?�t�    ?��4    >��    >w�?    >�~    >lL    Cb��    @U��    Cf/\    ?�\    A���    AP�)    @�     @�         C��    B��5    Cš�    C�NV    ?���    =�|    >X��    A�F    A��+    A�9�    A��h    A�|�    A�D    ?�	B    >ݸ    >���    ?�Z�    >��'    >��L    ?�I�    ?�v    ?���    �C(    ��t�    ���@    >^g�    C�\J    ?ɟ    C���    @P    A�S[    AQ+    @�}     @�}         C�Yy    B��x    C�iX    C�_    ���    ��E    �2�    A� �    A��l    A�@O    A���    A���    A��w    ��p�    �	�C    ��'�    ��>�    �
�*    ��(�    ?�T�    ?�x    ?��e    �B�    <��9    ����    >P1�    C���    ��8\    C��    @.|    AѾ    AO��    @�^     @�^         D"�    C T9    C�h�    C��    ��    ��cs    ����    A�'�    A��7    A�F�    A�ԕ    A�,�    A���    ��R    ��8�    ��P    ���    ��6�    ���    ?�[    ?�w�    ?��K    ���    �!c�    ���    >J�>    C�D{    ����    C���    @:�_    A�33    AO��    @�?     @�?         D)��    C �}    CΞ�    C��)    ���    �'_p    �$�    A�0�    A��7    A�M    A�J�    A�ѷ    A��N    ��Y�    �I    �ܛ(    ��d    �"    �ܜ�    ?�b�    ?�w�    ?���    �\]    ���    ���    >A��    C�B-    ����    C�P�    @;��    A즵    AP��    @�      @�          D>�Z    C	��    CϘ1    C�h    ��=    �c��    �,      A�<6    A��7    A�Q�    A��    A��    A��;    ���    �R��    ���    �˧H    �S��    ��M    ?�m]    ?�x    ?���    �3?�    ���     ��J    >8�    D�#    �+u    D	8!    @;�!    A    AP�"    @�     @�         DU6    C�    C�'�    C���    �Z�    ���l    �[�H    A�@�    A��Y    A�U�    A��p    A��^    A�Q�    �|P    ��2�    �=�    ���    ���"    �=�    ?�p�    ?�ud    ?��g    ���    ����    ��N'    >1�    Ds    �a��    D�-    @<�k    A���    AQV�    @��     @��         DaCT    C4�    C�pb    C�/;    �BY    ���E    ��I�    A�I�    A��o    A�]�    A��"    A��    A�L�    �-_p    ���    �lff    �.��    ����    �li�    ?�y)    ?�p�    ?���    ���@    ���    ��    >+��    D(
    ����    D#��    @<��    B }<    AP��    @��     @��         DqAX    C�R    CՔ�    C��     �r�c    ���    ��c     A�S&    A�|�    A�c�    A�<�    A��    A��4    �U_    ��=q    �� �    �V�!    ��l"    ��    ?큮    ?�l�    ?��.    ���g    ��G    � �    >(j�    D5Y�    ��>�    D0    @<�    Bބ    AP�'    @�     @�         Dd<    C8    C��%    C��    ����    ���M    ���    A�\�    A�tT    A�k�    A���    A�W
    A�:    �{�    �    ��8�    �}u�    ���    ��<    ?�3    ?�eV    ?�    ��     �ْ�    ��X    >&G4    D,5�    ����    D&�    @=�    B7�    AO��    @�B�    @�B�        Dr6V    C &�    C�`!    D �#    ���I    �	�    ���p    A�b    A�l"    A�r|    A�_�    A�2-    A��    ��%�    � �    ��    ��?�    � �N    ��S    ?��    ?�^     ?�J    ���p    ���    �&I�    >)>�    D6�X    ��j�    D/�y    @=�    B�D    AO��    @�     @�         D[��    CO�    C�$�    D ��    ��C�    �ɠ    ����    A�dZ    A�d&    A�v+    A�'�    A��    A���    ���    �Kx    ��4    ��:^    �6    ��8�    ?�    ?�V�    ?��    ���_    �׍P    �$�    >)7�    D%x�    ����    D`s    @=2    B
    AP��    @�#�    @�#�        DD	�    C    C�ļ    C�u�    ��    ���    ��    A�^�    A�Z    A�rG    A�x    A��    A�g    ���k    ��b    ��ӏ    ���    ���    ���    ?��    ?�N    ?��    ��L    ��%    �TL    >)�H    DOm    ��p�    DH�    @=�    B�    AP�6    @�     @�         DKe�    C˅    C�y7    D �/    ��O    �	>]    ����    A�RT    A�Ov    A�hs    A��}    A��2    A�f�    ��nc    ���    ���    ���!    �c�    ����    ?퀲    ?�D�    ?�q    ���    ��X�    ���    >(k�    D�^    ��t�    D*�    @=�    B
�    AO�v    @��    @��        D#<    B�5�    Cآ    C���    ����    ����    ���~    A�GE    A�F?    A�]d    A��V    A�r�    A�/O    ����    ���    ���P    ���K    ����    ����    ?�w2    ?�<`    ?��H    ����    ��z%    ���:    >.6�    C�ۦ    �ː�    C��    @=M    B�    AO��    @�u     @�u         D,;u    B�[#    C؍q    C��    ����    ��    ���f    A�1�    A�<6    A�K�    A�S    A��B    A�r�    ���/    ��    ��Vm    ��~    ��?H    ��Ws    ?�c�    ?�33    ?��    �ɤj    ��<�    ����    >1*    C�{�    ��e�    C���    @=*0    B;�    AP�    @��    @��        D$K�    B�9�    Cي�    C��    ��yr    ��k�    ���;    A��    A�1�    A�:*    A�gm    A��#    A�b    ���=    �Ƹ�    ����    ���0    ��ѷ    ����    ?�Q/    ?�*    ?��'    ��ӄ    �S;u    �(�    >8�    C�    ��wf    C�h    @=&B    BP}    AP:�    @�V     @�V         D��    B�7L    CېB    C���    ���+    ���    ��I�    A�    A�'�    A�,    A���    A���    A��    ����    ����    �~��    ���V    ����    �~��    ?�9�    ?� �    ?��`    ��oT    ����    �6�    ><|    C�o    ��f�    C�ɺ    @=F    Bo    AQ�    @�ƀ    @�ƀ        C��    B�)�    C��    C��u    ��M6    ��@�    �]YK    A��Q    A�!    A�C    A��    A�l�    A��	    ����    ��b�    �]�z    ���0    ��H�    �]��    ?�%F    ?�P    ?ٿH    �2C    �_�    <�QK    >D�6    C�r�    ���B    C��y    @<��    B^    AP�E    @�7     @�7         CӘ�    B���    C�Ͼ    C�
�    �p!    ��Y    �+�    A���    A��    A�	    A�S�    A��a    A���    �b��    ����    �(�    �d��    ��\�    �(      ?�p    ?�    ?ٮ>    �6�    ��;C    �A �    >L4;    C�NV    �e�|    C�S�    @<��    B F�    AP�E    @���    @���        C�`!    B�{q    C��    C�iX    �K~�    �]�    ���N    A���    A��    A��`    A���    A��    A�\]    �D0U    �\�1    �H�    �E�    �]��    �E9    ?��#    ?��    ?ٝ�    ��2    9Q��    ?%�?    >W��    C>��    �7    C2�    @<a    A��b    AP��    