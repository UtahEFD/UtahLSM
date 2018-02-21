CDF   0   
      time             command_line      sebs_ingest -s sgp -f E14      process_version       ingest-sebs-1.5-0.el6      ingest_software       ingest-sebs-1.5-0.el6      dod_version       sebs-b1-1.4    site_id       sgp    facility_id       E14: Lamont, Oklahoma      
data_level        b1     input_source      ?/data/collection/sgp/sgpsebsE14.00/Table30.20170626_000000.raw     resolution_description       The resolution field attributes refer to the number of significant digits relative to the decimal point that should be used in calculations. Using fewer digits might result in greater uncertainty. Using a larger number of digits should have no effect and thus is unnecessary. However, analyses based on differences in values with a larger number of significant digits than indicated could lead to erroneous results or misleading scientific conclusions.

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
datastream        sgpsebsE14.b1      history       Zcreated by user dsmgr on machine ruby at 2017-06-26 01:21:01, using ingest-sebs-1.5-0.el6         G   	base_time                string        2017-06-26 00:00:00 0:00       	long_name         Base time in Epoch     units         $seconds since 1970-1-1 0:00:00 0:00         Pt   time_offset                 	long_name         Time offset from base_time     units         'seconds since 2017-06-26 00:00:00 0:00          P�   time                	long_name         Time offset from midnight      units         'seconds since 2017-06-26 00:00:00 0:00          P�   qc_time                 	long_name         :Quality check results on field: Time offset from midnight      units         	unitless       description       vThis field contains bit packed values which should be interpreted as listed. No bits set (zero) represents good data.      bit_1_description         9Delta time between current and previous samples is zero.       bit_1_assessment      Indeterminate      bit_2_description         fDelta time between current and previous samples is less than the delta_t_lower_limit field attribute.      bit_2_assessment      Indeterminate      bit_3_description         iDelta time between current and previous samples is greater than the delta_t_upper_limit field attribute.       bit_3_assessment      Indeterminate      delta_t_lower_limit       @�        delta_t_upper_limit       @�,        prior_sample_flag               comment       �If the 'prior_sample_flag' is set the first sample time from a new raw file will be compared against the time just previous to it in the stored data. If it is not set the qc_time value for the first sample will be set to 0.         P�   down_short_hemisp                   	long_name         -Downwelling shortwave hemispheric irradiance       units         W/m^2      	valid_min                	valid_max         D�     
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
   missing_value         �<         Q�   qc_battery_voltage                  	long_name         0Quality check results on field: Battery voltage    units         	unitless       description       7See global attributes for individual bit descriptions.          Q�   lat              	long_name         North latitude     units         	degree_N       	valid_min         ´     	valid_max         B�          Px   lon              	long_name         East longitude     units         	degree_E       	valid_min         �4     	valid_max         C4          P|   alt              	long_name         Altitude above mean sea level      units         m           P�YPN�Bm�����C��                     BC��    A�v    C�)    C�o\    �8��    ��w2    ��[�    A���    A�u�    A�c    A�9X    A�8�    A�v�    �N	�    ����    ���    �P2�    ��D(    �	�!    ?��    ?�-�    ?�a    >�    ?]��    ?^i'    >7�X    ��    ��I    �'Y6    @8�>    A�nc    AQ��    @�      @�          C�    B��    Cϻ#    C뭑    �4    �H��    �h��    A�ɺ    A�v+    A�v�    A��.    A�e    A�ϫ    ���    �s�2    ��$_    �6�    �u2@    ��|�    ?��-    ?�.    ?��    �Z\    >2��    ?2@�    >b�    Bm��    ����    Bit�    @;T�    A��}    AQ˒    @�      @�          BӮ    A��~    C���    C�H    ��{�    �imn    �f�D    A���    A�y>    A�l�    A���    A�ɺ    A�H    ��͊    �c�    ��R    ��D�    �eY�    ��@d    ?��    ?�1    ?�	�    >�$�    ��xK    ?�s    >b�    A�ff    ���    A�)*    @6|�    A�oi    AQ�	    @�     @�         Bc��    Ah0U    CĨR    C�&F    ��RT    ?    �[��    A��?    A�|P    A�uZ    A��v    A���    A�	l    ��a�    ?�S�    �z��    ����    ?Ȕp    �|k�    ?��:    ?�3�    ?�4    �O��    >��    >�b    >�Cx    �HB�    ����    �I;d    @�Q    A��E    AQ�    @�      @�          @��    ?���    C��{    C�G
    �r��    @���    ?�:    A��'    A���    A�z�    A��    A�#:    A�F�    ����    @�C    ?��    ��_F    @�:    ?��    ?�֌    ?�8�    ?�    >��+    ?$�%    ?~˸    >C�    �U�g    ?�@:    �N�W    ?���    A�bN    AQ�    @��     @��         ��#O   �obc   C�D9    C�r�    >܌r    @�mH    @:��    A��    A���    A�uZ    A�[�    A��    A�X�    ���    @�;O    ?ܯ�    ����    @�K�    ?�#�    ?�ł    ?�<�    ?�4    ?�+    ?�    ?���    >V(�    ��    @C�:    �sd�    ?���    A���    AQ[W    @�     @�         AHɆ    �}W�   C��D    C�@�    @�    @Ѱ�    @r�f    A���    A���    A�p    A��q    A�Υ    A��    ?��A    @�a|    @��    ?�qL    @    @��    ?���    ?�@y    ?��    ?eM�    >�z�    ?��    ����   �<z    @�!    �+E9    ?�2    A��Q    AQy�    @Ȝ     @Ȝ         AD�"    ��gM   C��s    CИ�    @^    @�M    @��H    A��f    A���    A�a    A�L0    A�e�    A�'R    @�d    @�m�    @@�A    @,    @�҉    @Ał    ?��-    ?�At    ?��.    >�    ?Y�$    ?~�    ��p�   �2E�    @���    ���    ?���    A��o    AQc    @�      @�          A>8�    �`"6   C��H    Cψ�    @dq    A�    @��    A�~�    A���    A�M�    A���    A��    A�L�    @@�    @�jU    @h�    @B�f    A V    @i�?    ?���    ?�D    ?���    ?��    ?1��    ?Z��    ���;   �,"    @�ٔ    ��#n    ?�f�    A��	    AQ�q    @Ϥ     @Ϥ         A�    ��=�   C��
    C�ٺ    @���    AQn    @�y    A�r�    A���    A�AU    A��    A��    A�!�    @sHV    A��    @t��    @u�    A��    @v1'    ?��,    ?�F�    ?��     ?d�r    >�T    ?���    ��u�   �ҽ    @�dZ    ���    ?�j    A�;�    AQi�    @є     @є         ��P   ���<   C�H1    C�l)    @�c�    A�    @�\>    A�r�    A���    A�3h    A�T�    A�c     A��    @�B    @��!    @oc5    @���    A     @p��    ?��,    ?�H�    ?�֌    ?nO"    ?^ 6    ?��P    @ �   �N�    @�;O    �7�    ?�[�    A�u    AQ    @�V     @�V         A��q    �M�   C�O�    C��L    @�g�    A
r    @�S�    A�j    A���    A� '    A�Υ    A��H    A�P�    @�Ë    Ao    @��$    @�,(    A��    @�~�    ?���    ?�KI    ?�ł    >���    ?i�    ?^�"    ��[   ��q�    @�E�    ���D    ?�z�    A�^5    AQqv    @�     @�         B?��    ����   C��
    C˽�    @��    A�    @���    A�h>    A��=    A��    A�3�    A�+k    A��!    @�0�    A�    @���    @��y    A��    @���    ?���    ?�O7    ?�    ?:��    ?1�}    ?/7C    ���u   >��    @ՃQ    @╖    ?���    A� �    AP�K    @��     @��         A�|�    ��ߏ   C�,�    C�ɺ    @���    A�    @�i�    A�e,    A���    A��    A�|P    A��O    A���    @��A    A��    @�/Z    @�'�    A̎    @�%�    ?��    ?�R*    ?�C    ?;q�    ?1�@    ?�    ���   �X��    @�D=    ��M    ?���    A���    AP�    @؜     @؜         @�v    ����   C���    C�j    @��    A�+    @� �    A�]�    A��b    A��	    A��    A�x    A�+6    @�U    A\�    @��,    @��M    A	9�    @��M    ?��~    ?�T�    ?��    ?:�+    >���    ?c�    ���   ��    @ӭm    �c�    ?��x    A�>w    AQq�    @�^     @�^         A֡    �}��   C��    C�NV    @�P�    A��    @�ti    A�\)    A��:    A��    A�IR    A��    A��    @�e�    A�    @��V    @��    A�R    @��,    ?���    ?�Vm    ?雑    ?;
�    ?1��    ?^��    ��vS   ��Q    @���    ��R�    ?�|    A���    AO��    @�      @�          A�l�    �t�   C���    C�^�    @���    A�Y    @�?�    A�Ov    A��    A��2    A��6    A�7    A���    @�dE    A�o    @���    @���    Av�    @��g    ?��    ?�Z    ?鑼    ?0    ?V*    ?2�K    �NP�   �j     @�}�    @@��    ?�ٔ    A���    AQs�    @��     @��         B�޸    ��   C�f�    C�f    @�'�    A	�    @��`    A�PH    A���    A��    A�,q    A��    A�2a    @�{5    A(m    @��.    @�t    A I    @��X    ?��^    ?�[    ?郥    ?ل    >B�    >$    ���   B�\x    @Ĵ�    B̈́�    ?�     A�k�    AO�o    @ߤ     @ߤ         CQ    ��(�   C�-    C��y    @���    A4�    @�	    A�F�    A��    A��?    A���    A�5?    A��B    @��    AE$    @��A    @��n    A�    @�Z\    ?�w�    ?�^J    ?�u�    ?;k    ?1ދ    ?5u    ��X�   C2-    @�sm    C8�    ?�GZ    A��    AP��    @�     @�         B؆%    ���u   C���    C�u�    @�[B    A�    @�A_    A�F    A��    A���    A�9�    A��
    A�>�    @�ɛ    A��    @�a    @�X�    Ar2    @��H    ?�w2    ?�_�    ?�h�    ?]    ?5�u    ?b@�    �En�   B��b    @ϼ�    B�Y    ?�A�    A�u%    AQD�    @�     @�         A�&�    ����   C�\�    Cʙ    @���    A�r    @�.4    A�>    A���    A��q    A��=    A�Y�    A���    @�b�    AVM    @�,=    @���    A;�    @�'�    ?�pP    ?�c5    ?�^�    >���    �u�S    ��ƕ    �q��   �8�P    @Ș5    ����    ?�L�    A���    AP�K    @�u     @�u         B�6�    ���   C�iy    C�q    @���    A.�    @�v�    A�>�    A���    A���    A��    A�S    A�H�    @�sm    A
82    @��M    @���    A�    @��_    ?�p�    ?�d0    ?�Tv    ?j�    ?_g    ?�zN    �,�   B���    @���    B��    ?�9    A�_    APzx    @�V     @�V         B�1    ��g8   C�3    C�a    @���    A"Y    @�s�    A�>B    A��?    A���    A��k    A��x    A��t    @�G�    A�    @�\�    @��    A�j    @�[-    ?�pP    ?�f'    ?�M@    ?�    ?1�V    �\�    �LF*   B=o�    @��    BY�g    ?�e,    A�]/    AQ.I    @�7     @�7         @�R�    >��J    C�    C��#    @��    A4>B    @���    A�>�    A���    A���    A�C    A�    A�A�    @���    A(o    @���    @�n�    A)#�    @���    ?�p�    ?�hs    ?�<6    ?(    ?1�;    ?�R    =$�!    ���    @�ϖ    ��o�    ?��    A�>�    AP�$    @�     @�         B*��    @���    C�5    C�?�    @ز�    A8��    @�4�    A�9$    A��X    A�{    A���    A�_;    A�}V    @�h�    A,�    @�jj    @�b�    A-��    @���    ?�k�    ?�j    ?�3	    ?:�Z    ?1�Q    ?    >��    ����    Am�    >!^    ?�b�    A���    AP=�    @��     @��         C1    A�7    C�_;    C�T    @�@    Ao�    @�5�    A�,�    A���    A�s�    A��    A���    A��o    @�E9    Ap�    @�l�    @�3�    Ao�    @���    ?�`�    ?�la    ?�,|    ?>��    8�""    ?e�    >��    B�خ    @�(    B���    ?��    A�+k    AQ/    @��     @��         C�v�    B #T    C�	X    C�]�    @�Κ    A	U    @�^     A�.    A��    A�l"    A��_    A���    A�{    @�     A3�    @�ߏ    @��    A    @��d    ?�a�    ?�m]    ?�%F    ?Gk    >�{    <��w    =�n�    CP1'    @��~    CV��    ?�4�    A�>    AQ �    @�     @�         �<    Bf�    C��D    C�R    @�I�    @ϲ�    @k+    A�-    A���    A�`B    A�+    A�Vm    A�q�    @���    @�d    @r�    @��    @β�    @t0�    ?�`�    ?�p�    ?��    �9�    <���    �]    �<    �<    @���    �<    ?�i    A��    AQ(    @�     @�         �<    B���    C��s    C�5�    @���    @���    @)�H    A�2a    A�Ĝ    A�O�    A��%    A�J�    A�t    @s�    @�G�    @;�	    @u��    @�D�    @<�5    ?�e�    ?�s�    ?�    >|��    >q�    ��~�    �<    �<    @|�    �<    ?�N�    A��5    AQ��    @�}     @�}         �<    B��!    C�N5    C�p�    @:C-    @/�~    ?��]    A�1�    A��-    A�(�    A��B    A���    A��    @4&    @+��    ?�/    @6�    @,�>    ?��    ?�e�    ?�r�    ?��O    =�)W    ==�S    ��>�    �<    �<    @��    �<    @H    A�:^    AO�s    @�^     @�^         Dg�    B�}�    C��Z    C�xs    ?��    �T�&    �ݓ�    A��    A��H    A�
	    A���    A���    A�n/    ?ӑ�    <�S    >e �    ?տH    <�$    >f�]    ?�O�    ?�oT    ?��    ���7    �Y�J    �(l�    >Bӄ    C��B    <�=�    C���    @/]    A�!�    AQ�[    @�?     @�?         D!�    B��T    C�|)    C�0�    ?6�z    ���    � ��    A�
=    A���    A��    A���    A�P    A�    >��    ���U    ��^�    >���    ��|�    ��k�    ?�B    ?�g�    ?行    >nfE    =m    �%7�    ><��    C�u`    ���    C���    @4��    A�~(    AP�    @�      @�          DH"N    C��    C�p     C�z    ��	W    ��>�    ���Y    A��    A���    A��b    A�`B    A�J�    A���    �U9"    ���    �^�/    �Wg^    ���    �_��    ?��    ?�\�    ?�q    �4�q    �-&�    ����    >8L    D	�'    �tC    D��    @8�>    A��    AO��    @�     @�         DX_m    C3    Cž�    C���    �=P]    ���    ��¤    A��R    A���    A�rG    A��    A��]    A�I    ���    ��o*    ���P    �8    ���     ��2a    ?�    ?�P�    ?�G    �5 �    � �l    ��@�    >5i    D�    ��}A    Dw    @;    A���    AQ:�    @��     @��         DTl�    C�    C�N�    C��    ��]�    �44    �+kQ    A��    A��    A�P    A�rG    A��i    A���    �u��    �(R�    �
�    �x    �)bN    ��    ?�}    ?�H    ?�(�    �
��    �-�    ��Ft    >4��    D}�    ���    DK�    @;Y`    A�    AQ2�    @��     @��         DM��    C��    C̝    DD9    ���m    �qG�    �P/    A��)    A��~    A�#�    A�<�    A�/    A�C�    ���{    �_	    �6�I    ����    �`p;    �7�    ?�
    ?�B    ?�P    ���    ����    ��5    >68�    DU    �9��    Dr�    @;��    A��a    AP�    @�     @�         D�L�    CDr�    C���    D�D    ��^�    ��po    ��S    A��6    A��    A���    A�*e    A�3�    A��    ��$�    ���    �_�m    ��nn    ����    �a�    ?��    ?�:�    ?���    �_��    ���    ��C    >;=    D=-P    �d�h    D9��    @=    B��    AQ2a    @�B�    @�B�        Dx��    C2�+    C��7    DC�    �*�    ��N    ���    A���    A�~(    A��    A�*�    A�a|    A��[    ���    ����    ����    ��    ���    ����    ?�!�    ?�5T    ?��    ����    ��ϖ    ���o    >7�:    D,>5    ��x    D'��    @=-#    Bu    AQ_p    @�     @�         DVAh    C�w    C�|j    D��    �B�=    ��c�    ���a    A��    A�q�    A���    A�r    A��z    A�    �)�t    ����    ���0    �+q�    ����    ����    ?�d    ?�*�    ?�}    ��J�    ��\    ���    >7��    D��    ��!-    D�    @<�b    B�    AQ�    @�#�    @�#�        DPKT    C��    C��    D/\    �KIR    ��.    ��JX    A��    A�g8    A���    A���    A�(�    A�y>    �2A�    ��k�    ���k    �4�    ���h    ��r    ?�    ?�!    ?�s    ���    ���    ����    >:N{    D
Bo    ���v    D��    @=K    B��    AP�6    @�     @�         DP�    Ct�    C�)y    D*o    �Y
�    ��+�    ��|    A��B    A�`v    A�b    A�
r    A�g8    A�    �?�6    ��J�    ���    �A��    ��~]    ����    ?��    ?�    ?�T�    ��Jw    ��پ    ���e    >@	p    D�L    ��dZ    D�J    @=W    B	��    AR@�    @��    @��        DP�s    C#��    Cͺ    D	HB    �qB�    ��Dg    ��	l    A��    A�X    A�C�    A�lW    A���    A�r�    �UA     �ʕ�    ���&    �Wo�    ��خ    ���,    ?�     ?��    ?�:    �Η9    ����    ��TL    >H��    D��    ���D    D (�    @=h�    Bf�    AP��    @�u     @�u         DJ�R    C)-    CЈ�    D	ى    �~;�    ����    ���    A�՛    A�O    A��    A��    A���    A�H�    �^��    ��E9    ��:�    �aC�    �ɄM    ��2�    ?�@    ?�c    ?�b    ���    ��P    ��    >U�I    C�{#    ��2�    C�    @=rG    BV    AO�s    @��    @��        D3)    CnV    C��    D	d�    �v�M    ���    ���'    A��K    A�G�    A��v    A�q�    A��    A�'�    �Z��    ���(    ����    �\��    ��*�    ����    ?�t    ?�)    ?��    ��A�    ���K    ��O    >^J8    Cڪ    ���    C��=    @=fQ    B��    AQ�a    @�V     @�V         C��h    Bu�    C��    D GL    �j�    ����    ����    A���    A�@�    A���    A��d    A�0�    A�Υ    �T/�    ���#    ��6F    �VX�    ��ں    ���    ?��U    ?���    ?��9    ���}    ��V    ��E�    >S�>    C��    ��&    C�    @<@�    B��    AQ?    @�ƀ    @�ƀ        Cp�d    BD��    C׈1    C�[D    �A��    ����    �_��    A���    A�<6    A���    A��    A��E    Aƴn    �1�    �qJ�    �[w2    �2��    �rɆ    �\��    ?�ٔ    ?���    ?�    �j�    �i��    �@
    >Q	�    B�Z�    �`�    B��}    @:��    B �M    APM�    @�7     @�7         C�    B]AU    C�
=    C�H�    �$�F    �F�e    �%��    A�rG    A�:�    A��    A���    A�N�    A��;    ��v    �B�F    �'��    �N'    �C�    �(e,    ?���    ?���    ?楤    �t_[    �/�R    >G�    >]��    C�    �0bN    B��Z    @9s�    A��S    AQ�    @���    @���        CMH�    B12�    C��H    C�o�    ���    �/��    ��    A�X�    A�9�    A��M    A�F?    A�    AƜ    ��_�    �)�    ���    � xb    �*+�    �	�J    ?��+    ?��b    ?��    ���6    ����    >1'    >]T    B�F�    �&l    B�!�    @9�    B W?    AR�    