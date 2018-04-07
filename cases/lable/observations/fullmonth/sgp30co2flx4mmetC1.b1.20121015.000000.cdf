CDF   0   
      time             site_id       SGP    facility_id       C1     sds_mode      b1     missing_value         -9999      
sample_int        0.12 seconds   averaging_int         
30 minutes     instruments       0Anemometer:Fill Windmaster Pro, IRGA: LICOR 7500   time_stamp_details        RThe time stamp currently reflects the _start_ of the 30 minute sampling interval.      boom_direction        K0 degrees; wdir is accurate to within an offset of approximately 10 degrees    sign_convention      fc_corr, h, le, mean_g_soil are positive away from surface, r_net is positive towards the surface. That is r_net = h+le+g. Also, r_net = (r_down_short_hemisp+r_down_long_hemisp) - (r_up_short_hemisp+r_up_long_hemisp). fc_corr is defind such that (fc_corr-fc_wpl_h-fc_wpl_le) = fc_uncorr     contact_name      %Marc Fischer, LBNL, mlfischer@lbl.gov      PLEASE_READ_THIS     �For all publications and presentations, please acknowledge: 'U.S. Department of Energy as part of the Atmospheric Radiation Measurement Program.' The automatic inclusion of a data originator as a co-author is not insisted upon in the ARM program, but the source of any data should be clearly recognized either as a co-author or through an appropriate acknowledgment.  The ARM/LBNL Carbon Project contact for this data set is Marc Fischer (mlfischer@lbl.gov). Also please note that we may make adjustments to the data to incorporate adjustments to calibration scales or other issues.  Users should contact Marc Fischer (mlfischer@lbl.gov) to inquire about planned data releases and revisions. Whenever possible, we would appreciate receiving preprints for publications that use the data to insure that the quality and limitations of the data are accurately represented. Your questions and comments are welcome.      MENTOR_QC_FIELD_INFORMATION      �The qc flags use binary, or bitmapped, values (ie values of 0,1,2,4,8,etc.) to note each relevant condition. These can be summed for more detailed qc information. The qc values are:
QC VALUE:  EXPLANATION:
 0-         value not suspect
 1-         value missing
 2-         value below minimum or above maximum or value is +-infinity
 4-         one or more dependencies failed; see the 'dependency' attribute of variable; if a dependency is a 
            spike count, it fails if it is > 100, if a dependency is a variable, it fails when the variable's qc flag
            is not equal to 0
 8-         value has large variance; see 'large variance condition' attribute of variable
16-         value susbect because it has more than 100 spikes; see 'dependency' attribute of variable;
            currently this flag applies to t,q,c
32-         value suspect because of another condition, see 'special condition' attribute of variable.
            Currently only applies to fc_corr, ustar.   	arm_field               irga_serial       C     sonic_serial      B     nrlite_serial         B     history       �created by code 4msonicb1tob1met.c, version r11, with operating system RedHat Linux, kernel 2.4.18-18.7.x, i686 on Apr 09 2015, 08:44:30 GMT      k   	base_time                string        15-Oct-2012,00:00:00 GMT       	long_name         Base time in Epoch     units         #seconds since 1970-1-1 0:00:00 0:00         M�   time_offset                 	long_name         Time offset from base_time     units         'seconds since 15-Oct-2012,00:00:00 GMT          M�   yyyydddhhmmss                   units         yyyydddhhmmss      	long_name         start of integration interval           M�   doy                 units         fractional days    	long_name         fractional day of the year          M�   fc_corr                 units         umol m-2 s-1   	long_name         WPL corrected CO2 flux     	valid_min         ����   	valid_max               
dependency        !nspk_unrot_w,mean_t,mean_q,mean_c           M�   
qc_fc_corr                  units         unitless   	long_name         qc flag for fc_corr         M�   fc_wpl_h                units         umol m-2 s-1   	long_name         #additive WPL H correction to c flux    	valid_min         ����   	valid_max                    M�   qc_fc_wpl_h                 units         unitless   	long_name         qc flag for fc_wpl_h        M�   	fc_wpl_le                   units         umol m-2 s-1   	long_name         $additive WPL LE correction to c flux   	valid_min         ����   	valid_max                    N    qc_fc_wpl_le                units         umol m-2 s-1   	long_name         $additive WPL LE correction to c flux        N   h                   units         W m-2      	long_name         corrected sensible heat flux   	valid_min         ���8   	valid_max               
dependency        nspk_unrot_w,mean_t         N   qc_h                units         unitless   	long_name         qc flag for h           N   le                  units         W m-2      	long_name         WPL corrected latent heat flux     	valid_min         ����   	valid_max               
dependency        nspk_unrot_w,mean_t         N   qc_le                   units         unitless   	long_name         qc flag for le          N   
mean_rot_u                  units         m s-1      	long_name         mean horizontal wind speed          N   qc_mean_rot_u                   units         unitless   	long_name         qc flag for mean_r          N   mean_t                  units         degree C   	long_name         5mean sonic temperature (t), i.e. virtual temperature       	valid_min         ����   	valid_max            2   large_variance_condition      abs(mean_t)/sqrt(var_t) < 2         N    	qc_mean_t                   units         unitless   	long_name         qc flag for mean_t          N$   mean_q                  units         mmol m-3   	long_name         mean H2O density   	valid_min         <#�
   	valid_max         D���   large_variance_condition      abs(mean_q)/sqrt(var_q) < 2         N(   	qc_mean_q                   units         unitless   	long_name         qc flag for mean_q          N,   mean_c                  units         mmol m-3   	long_name         mean CO2 density   	valid_min         A (�   	valid_max         A��        N0   	qc_mean_c                   units         unitless   	long_name         qc flag for mean_c          N4   mean_p                  units         kPa    	long_name         mean IRGA pressure     	valid_min            ^   	valid_max            e        N8   	qc_mean_p                   units         unitless   	long_name         qc flag for mean_p          N<   	var_rot_u                   units         (m s-1)2   	long_name         variance of u      	valid_min                	valid_max            d        N@   qc_var_rot_u                units         unitless   	long_name         qc flag for var_r           ND   	var_rot_v                   units         (m s-1)2   	long_name         variance of v      	valid_min                	valid_max            d        NH   qc_var_rot_v                units         unitless   	long_name         qc flag for var_r           NL   	var_rot_w                   units         (m s-1)2   	long_name         variance of w      	valid_min                	valid_max                    NP   qc_var_rot_w                units         unitless   	long_name         qc flag for var_r           NT   var_t                   units         C2     	long_name         variance of t      	valid_min                	valid_max                    NX   qc_var_t                units         unitless   	long_name         qc flag for var_t           N\   var_q                   units         (mmol m-3)2    	long_name         variance of q      	valid_min                	valid_max           p        N`   qc_var_q                units         unitless   	long_name         qc flag for var_q           Nd   var_c                   units         (umol m-3)2    	long_name         variance of c      	valid_min                	valid_max         >���        Nh   qc_var_c                units         unitless   	long_name         qc flag for var_c           Nl   wdir                units         degrees    	long_name         horizontal wind direction      	valid_min                	valid_max           h        Np   qc_wdir                 units         unitless   	long_name         qc flag for wdir        Nt   theta                   units         degrees    	long_name         rotation to zero w     	valid_min                	valid_max           h        Nx   qc_theta                units         unitless   	long_name         qc flag for theta           N|   phi                 units         degrees    	long_name         rotation to zero v     	valid_min                	valid_max           h        N�   qc_phi                  units         unitless   	long_name         qc flag for phi         N�   ustar                   units         m s-1      	long_name         friction velocity      	valid_min                	valid_max            
        N�   qc_ustar                units         unitless   	long_name         qc flag for ustar           N�   Lmoni                   units         meters     	long_name         Monin-Obukhov length scale     	valid_min         ����   	valid_max           #(        N�   qc_Lmoni                units         unitless   	long_name         qc flag for Lmoni           N�   bar_pres                units         kPa    	long_name         barometric pressure    	valid_min            ^   	valid_max            e        N�   qc_bar_pres                 units         unitless   	long_name         qc flag for bar_pres        N�   t_air_upper                 units         degree C   	long_name         !mean air temperature, upper level      	valid_min         ����   	valid_max            2        N�   qc_t_air_upper                  units         unitless   	long_name         qc flag for t_air_upper         N�   t_air_lower                 units         degree C   	long_name         !mean air temperature, lower level      	valid_min         ����   	valid_max            2        N�   qc_t_air_lower                  units         unitless   	long_name         qc flag for t_air_lower         N�   rh_upper                units         %      	long_name         #mean relative humidity, upper level    	valid_min                	valid_max            d        N�   qc_rh_upper                 units         unitless   	long_name         qc flag for rh_upper        N�   rh_lower                units         %      	long_name         #mean relative humidity, lower level    	valid_min                	valid_max            d        N�   qc_rh_lower                 units         unitless   	long_name         qc flag for rh_lower            N�   z_upper              units         meters     	long_name         height of upper RH/T sensor         M�   z_lower              units         meters     	long_name         height of lower RH/T sensor         M�   mean_g_soil                 units         W m-2      	long_name         -average of up to four soil heat flux sensors       	valid_min         ����   	valid_max            �   special_condition         3less than two out of four sensors read valid values         N�   qc_mean_g_soil                  units         unitless   	long_name         qc flag on mean_g_soil          N�   stderr_mean_g_soil                  units         none   	long_name         standard error of mean_g_soil           N�   z_g_soil             units         meters     	long_name         depth of soil heat flux probes          M�   r_net                   units         W m-2      	long_name         net solar radiation (NR-lite)      	valid_min         ����   	valid_max                    N�   qc_r_net                units         unitless   	long_name         qc flag for r_net           N�   r_tot                   units         W m-2      	long_name         5total downwelling shortwave radiation (LI-200 sensor)      	valid_min         ����   	valid_max           �        N�   qc_r_tot                units         unitless   	long_name         qc flag for r_tot           N�   r_down_short_hemisp                 units         W m-2      	long_name         -downwelling shortwave hemispheric irradiance       	valid_min         ����   	valid_max           �        N�   qc_r_down_short_hemisp                  units         unitless   	long_name         qc flag for r_down_short_hemisp         N�   r_up_short_hemisp                   units         W m-2      	long_name         *upwelling shortwave hemispheric irradiance     	valid_min         ����   	valid_max           �        N�   qc_r_up_short_hemisp                units         unitless   	long_name         qc flag for r_up_short_hemisp           N�   r_down_long_hemisp                  units         W m-2      	long_name         <downwelling longwave hemispheric irradiance; t_rad corrected   	valid_min         ����   	valid_max               
dependency        t_rad           N�   qc_r_down_long_hemisp                   units         unitless   	long_name         qc flag for r_down_long_hemisp          N�   r_up_long_hemisp                units         W m-2      	long_name         :upwelling longwave hemispheric irradiance; t_rad corrected     	valid_min         ����   	valid_max               
dependency        t_rad           N�   qc_r_up_long_hemisp                 units         unitless   	long_name         qc flag for r_up_long_hemisp        N�   t_rad                   units         degree K   	long_name         5radiometer body temperature (for longwave correction)      	valid_min            �   	valid_max           J   
dependency        t_rad           N�   qc_t_rad                units         unitless   	long_name         qc flag for t_rad           O    ppfd                units         umol m-2 s-1   	long_name         (photosynthetic photon flux density (PAR)   	valid_min         ����   	valid_max           �        O   qc_ppfd                 units         unitless   	long_name         qc flag for ppfd        O   mean_t_soil_upper                   units         degree C   	long_name         -average of two upper soil temperature sensors      	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O   qc_mean_t_soil_upper                units         unitless   	long_name         qc flag for mean_t_soil_upper           O   stderr_mean_t_soil_upper                units         none   	long_name         Sstandard error of mean_t_soil_upper; note maximum of two values used in calculation         O   mean_t_soil_middle                  units         degree C   	long_name         .average of two middle soil temperature sensors     	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O   qc_mean_t_soil_middle                   units         unitless   	long_name         qc flag for mean_t_soil_middle          O   stderr_mean_t_soil_middle                   units         none   	long_name         Tstandard error of mean_t_soil_middle; note maximum of two values used in calculation        O    mean_t_soil_lower                   units         degree C   	long_name         -average of two lower soil temperature sensors      	valid_min         ����   	valid_max            2   special_condition         (not one of two sensors read valid values        O$   qc_mean_t_soil_lower                units         unitless   	long_name         qc flag for mean_t_soil_lower           O(   stderr_mean_t_soil_lower                units         none   	long_name         Sstandard error of mean_t_soil_lower; note maximum of two values used in calculation         O,   mean_m_soil_upper                   units         cm3 cm-3   	long_name         +average of four upper soil moisture sensors    	valid_min                	valid_max         ?      special_condition         /less than two of four sensors read valid values         O0   qc_mean_m_soil_upper                units         unitless   	long_name         qc flag for mean_m_soil_upper           O4   stderr_mean_m_soil_upper                units         none   	long_name         $standard error of mean_m_soil_upper         O8   mean_m_soil_lower                   units         cm3 cm-3   	long_name         +average of four lower soil moisture sensors    	valid_min                	valid_max         ?      special_condition         /less than two of four sensors read valid values         O<   qc_mean_m_soil_lower                units         unitless   	long_name         qc flag for mean_m_soil_lower           O@   stderr_mean_m_soil_lower                units         none   	long_name         $standard error of mean_m_soil_lower         OD   precip                  units         
mm 30min-1     	long_name         mean precipitation     	valid_min            ^   	valid_max            e        OH   	qc_precip                   units         unitless   	long_name         qc flag for precip          OL   lat              units         degrees    	long_name         latitude of instrument location         M�   lon              units         degrees    	long_name          longitude of instrument location        M�   alt              units         meters     	long_name         #height of tower base from sea level         M�   ppfd_up                 units         umol m-2 s-1   	long_name         2Upwelling photosynthetic photon flux density (PAR)     	valid_min         ����   	valid_max           �        OP   
qc_ppfd_up                  units         unitless   	long_name         qc flag for ppfd_up         OT   zm               units         	meters         	long_name         #height of instrument from tower bas         M�   zrad             units         meters     	long_name         height of radiation instruments         M�   z_t_soil_upper               units         meters     	long_name         depth of upper t soil probes        M�   z_t_soil_middle              units         meters     	long_name         depth of middle t soil probes           M�   z_t_soil_lower               units         meters     	long_name         depth of lower t soil probes        M�   z_m_soil_upper               units         meters     	long_name         depth of upper m soil probes        M�   z_m_soil_lower               units         meters     	long_name         depth of lower m soil probes        M�P{R @|(�@	���L��Bk����jC�� @���@   �L�;����  ���;L��        B}H\V$  C�� =ǭ%  ����  ���e  ��&�  ���h  �?�0A �A�y2  �CՁ6  �As)  �B��N  �>�m  �=.�!  �<��  �>�0  �@�J�  �;�.�  �C�E
  �?���  �C$�  �=rĘ  �?ӹ�  �B���  �A<��  �A�33  �B9\)  �A��\  ���O�    @`�]  ��<  ���z�  �?Qhs  �C��?  �C��  �C��f  �?P �  �A���  �>�`A�=p  �>���A��\  �=aH <�E�  ��< �<  !��<       �      �@�      B}H\V߀ C���>;�W  ���^�  �:tg&  ��E  ���p-  �?X �A�m  �Cסo  �As�!  �B��$  �<�	�  �>��W  �:�K  �=��  �>Y#5  �:jvb  �C�+�  �?N��  �B̯-  �=E'�  �?��  �B��  �A.�H  �A�(�  �BH{  �Bff  �����    ?�p��J��  ��<  ����  �@"�!  �C��  �C�Pt  �C��3  �      �Al=p  �>޸PAj��  �?��Am�  �<u� <|z:  ��< �<  !��<       �      �@�      B}H\X�  C��U��BE  ��Ǯ  �����  ��gDw  ��:0^  �?ʜ�  �A�H�  �C�D:  �At�^  �B��c  �<��  �<�zO  �:z�d  �=a�@  �@MU�  �;#��  �Cz'�  �@�  �B�OM  �<��I  �?%��  �B��N  �A!��  �A���  �BT=q  �Bz�  ���t�    ?�lH�C�  ��<  ���y  �?l�D  �C�$�  �C�g%  �C��  �      �AW33  �>.�AV{  �>�\ AX�R  �=�
@<4cI  ��< �<  !��<       �      �@�     B}H\YP� C�� ��R^  ���  ��j�  ��;�A  ���r  �?�U  �A�Q  �C؆J  �Aw!  �B��t  �<�N�  �=�O�  �:�6  �>W�  �?���  �:���  �Ct>�  �?���  �B�}=  �;���  �<i��  �B��  �A
�R  �A��  �Ba33  �B%��  �����    ?|k��@�  ��<  ���P  �?�M�  �C�F-  �C�u  �C�s3  �      �AK�
  �>�AJff  �>u AM  �>8Q�;�e�  ��< �<  !��<       �      �@�      B}H\[  C������  ��7xD  �;�;�  ����+  �>��  �?���  �A�7  �C��P  �Av��  �B���  �=�Ѵ  �<�5�  �;�x�  �>9��  �@A�"  �;j'  �C`��  �?�E  �B2?�  �=��F  �C 94  �B��-  �A�
  �A���  �B`�  �B#
=  ����    ?w���>  ��<  ���+  �@|�  �C�*s  �C��  �C�ff  �      �A>Q�  �>   A=�  �>p��A@z  �=�
@;��8  ��< �<  !��<       �      �@��     B}H\[�� C��U�5�&  ��Ť�  �;z�@  ��
D  ��l�  �@(�b  �A�sa  �C�oU  �A|�  �B��N  �=�ŷ  �=!�+  �;�2%  �?j  �B$^�  �<��  �CV0W  �?���  �B�[  �=k�$  �@ ��  �B��v  �@��y  �Apz�  �B���  �B>
=  ���p�    ?ً1�L    ��<  ��7
=  �����  �C��;  �C��1  �C��  ��o  �A?G�  �>G� A>�\  �>k� A@Q�  �=L̀;��  ��< �<  !��<       ��D��  �@�     B}H\]w  C�� ����  ����!  �� �2  ���i  ��0o�  �@�  �A��a  �C��  �A$q  �B���  �=7�  �<�X�  �;��  �>&v  �@>�  �;]N�  �CV�M  �?�o#  �Bi1  �=V��  �?ąx  �B��  �@�G�  �Ak�  �B��R  �BE�H  ���    ?�9W�@33  ��<  ��͑h  �>l�D  �C�,�  �C�5�  �C��f  �      �A6��  �>G� A6    �>��`A7�  �=��;��N  ��< �<  !��<       ��D��  �@Ȝ     B}H\^2� C����)�W  ���
\  ��.y?  ����\  ���k  �@H�  �A�}j  �C�s  �A��?  �B���  �>*�(  �=q��  �;�V  �>�U�  �BHĈ  �<�u  �CPT�  �?�^�  �A�C  �=x�5  �@[�  �B��  �@��
  �Ab�H  �B�33  �BQ�  ����]    ?��Y�RQ�  ��<  ���\)  �?�1'  �C�Pt  �C�>�  �C�ff  ��o  �A>{  �=�A=�  �>�`A>��  �=8Q�;��  ��< �<  !��<       �      �@�      B}H\_�  C��U�\|C  ��!��  �;ßj  ��Z��  ��F  �@8/  �A���  �C�T  �A���  �B���  �>�  �=;��  �<z�+  �>���  �Ai  �<��  �CT.�  �?��  �B ��  �=��7  �@��\  �B�+  �@��;  �AUG�  �B�ff  �Ba=q  ���x�    ?��H�S  ��<  ��lI�  �?V  �C���  �C��  �C�&f  ��o  �AA��  �=�Q�AAp�  �=���AB(�  �<�� <�j  ��< �<  !��<       ��o  �@Ϥ     B}H\`�� C�� ��)  ��
�)  ��˨�  ��:0�  ���?Y  �@>��  �A�y  �C쏏  �A��  �B��  �=�M`  �<��f  �<�  �>4u�  �@�=8  �;��  �CR�n  �?��w  �A��r  �=�\  �@е  �B��q  �@�V  �AL��  �B�ff  �Bm��  ���
>    ?���P�  ��<  ��͑h  �?�l�  �C��  �C�|�  �C���  ��u  �A6=p  �=�Q�A6z  �=���A6��  �=\�<!  ��< �<  !��<       ����
  �@є     B}H\bY  C�����,~  �����  ��f�6  ���+�  ��S\�  �@:u  �A��Z  �C��n  �A���  �B��[  �=#;�  �=�<  �;ϕ�  �>i�  �A��5  �<W�  �CQ<P  �?��  �A��|  �=n�J  �@��  �B���  �@���  �AF�\  �B�    �Bw(�  ���Q�    ?���K(�  ��<  ���  �?��y  �C��  �C�*  �C��3  ���j  �A2Q�  �=�@A2(�  �=@A2�H  �<�� <�r  ��< �<  !��<       ��#�
  �@�V     B}H\c� C��U��-�  ��]&�  �8�T7  ���	�  ���||  �@?�  �A���  �C�X"  �A��  �B���  �=i�  �=��  �<d�  �>C�  �@�I�  �;���  �CQ�U  �?���  �A�
�  �=���  �@p\�  �B��'  �@n�R  �A@    �B�ff  �B{33  ���l�    ?���KQ�  ��<  ����R  �?�^  �C�?+  �C���  �C�s3  ��o  �A,Q�  �=��@A,(�  �=�R A,�H  �=\ ;��  ��< �<  !��<       ��t�  �@�     B}H\d�  C�� �S��  ��Cg�  ��f��  ���,  ����T  �@T٨  �A�p�  �C���  �A�M  �B���  �=��W  �=�R  �<��<  �=���  �@˿�  �;ZD  �CS��  �?��U  �A�ǚ  �=�7/  �@�̿  �B��`  �@[��  �A=�  �B���  �B�  ��� �    ?�IN�X=q  ��<  �����  �@�
  �C�	i  �C���  �C�s3  ��Y�  �A/�  �=�\@A/��  �=���A0=p  �<�� ;�L  ��< �<  !��<       ��D��  �@��     B}H\e�� C�����W�  ����8  ���8  ��܏(  ����  �@X��  �A��I  �C���  �A�zS  �B�ɿ  �>���  �=���  �=  �=���  �?�h  �;|1�  �CVQ�  �?��6  �B	F�  �>/�  �A3�l  �B��6  �@M��  �A9G�  �B���  �B�
=  ���+    ?�S��W=q  ��<  ���Ĝ  �@j  �C���  �C�E!  �C�@   ���o  �A.ff  �=�\@A.Q�  �=���A.�H  �<�� <VX  ��< �<  !��<       ��ix�  �@؜     B}H\g;  C��U�Ƣ�  ���Ǵ  ��>  �����  ����j  �@��C  �A���  �C��V  �A�fp  �B��  �>Vj�  �=�;�  �=Tv  �=� �  �?\6�  �;J��  �CT}  �?�)A  �B ]�  �>8�[  �A�[�  �B��)  �@a�  �A@��  �B�    �B~�  ��r^4    @ ���g=q  ��<  ���ȴ  �@)��  �C�_/  �C�<  �C���  ��H�9  �A6��  �=��A6��  �=��A7G�  �<u� <R^�  ��< �<  !��<       ��t�  �@�^     B}H\g�� C�� ��b�  ���i2  �:�k  ����?  ��Ê"  �@�ԏ  �A��{  �C�K�  �A�i
  �B��  �>�I�  �=��  �=AK  �=��k  �>��   �;j�  �CS�d  �?ѻi  �A��  �>3�c  �A��)  �B�ؓ  �@W;d  �A?�  �B���  �B  ����T    ?�KX�c33  ��<  ��
=  �?�?}  �C��  �C��  �C���  ��C�  �A6(�  �=aG�A6z  �=uA6z�  �<u� <u  ��< �<  !��<       �����  �@�      B}H\i�  C���>���  �����  �;��  ���X>  �����  �@rW%  �A�h�  �C�  �A�vh  �B��=  �>�  �=~�
  �<�͵  �=��f  �>��E  �:=x�  �CR��  �?�Z  �A��$  �>
��  �A-�  �B��7  �@H�  �A9�  �B���  �B��  ����l    ?�w��[��  ��<  ���  �?�ƨ  �C�x:  �C���  �C�ff  ��ix�  �A.z  �=uA.    �=��A.z�  �<u� <z  ��< �<  !��<       ����  �@��     B}H\jg� C��U>B��  ��d/F  �;�tV  ����  ����	  �@@��  �A�/�  �Cﶼ  �A�E�  �B��x  �>���  �=�b�  �<͆&  �=��q  �@H�  �:"c-  �CV_�  �?��  �B	~�  �=���  �@[i�  �Bý�  �@3"�  �A4    �B�33  �B�=q  ����8    ?��"�T��  ��<  ����  �?Ұ!  �C�pd  �C��  �C�33  ���C�  �A'�  �=���A'�  �=��A(z�  �<�� <|�  ��< �<  !��<       ��T��  �@ߤ     B}H\l  C�� �۫|  ���  ��/�O  ��Eʞ  ���a�  �?�	�  �A�R-  �C��  �A�J<  �Bþ�  �>f.&  �>M�  �<v��  �=��  �@֙  �:��  �Cb�f  �?�H  �B;ݕ  �=ȹ�  �@��R  �Büj  �@dZ  �A,    �B���  �B��
  ��5?    ?���7  ��<  ����  �?��y  �C�Vc  �C�6c  �C��   �����  �A�R  �>#� A�\  �>8R A�
  �<#� <UQ�  ��< �<  !��<       ���\)  �@�     B}H\l؀ C���>Lk*  ��^F�  ���X  ���Y�  ��E�  �?�<B  �A�ۦ  �C�yF  �A��  �Bû  �=��  �=d_%  �<^�  �=��  �@fq�  �:Y��  �C�d  �>�6�  �BՐ  �=<�  �?�Qv  �Bù�  �@.  �A4z�  �B���  �B���  ���M�    @6y�"\)  ��<  �      �@�|�  �C��X  �C���  �C��f  ��$�/  �A=p  �>��Ap�  �>���A��  �=\�<S�  ��< �<  !��<       ���h  �@�     B}H\n�  C��U��T]  ��=�  ���B�  ���X�  ��6�W  �@�-  �A���  �C��&  �A���  �Bý�  �=��_  �=ȷu  �<y-�  �=�.P  �@�#�  �;z�R  �C���  �?Ф�  �B�*$  �=��H  �@˒�  �Bú�  �@W�  �A@��  �B�33  �Bx  ��6E�    @4#�7=q  ��<  �      �@�+  �C��:  �C�[�  �C��  ��$�  �A%34  �>ff`A%33  �>\) A&    �<�� <h�  ��< �<  !��<       ��I�  �@�u     B}H\oI� C�� ��1_  �=~��  ����i  �?�/�  ����  �?�P>  �A�"�  �C��  �A�[�  �B��  �=\��  �>  �<�C  �>��  �?��  �;a �  �C���  �>��  �B�+W  �=�  ��W.{  �B���  �@�j  �ABff  �B���  �Bt�
  ��"z�    @��*�\  ��<  ����T  �@��  �C�G�  �C���  �C�ff  ���-  �A!��  �>aG�A!\)  �>ff�A"{  �<�� <p�  ��< �<  !��<       ���v�  �@�V     B}H\p�  C����	}>  ��g  ����  ��-y  ���Y[  �?���  �A��$  �C��  �A��1  �B��E  �<�h�  �>z C  �;���  �>Ȟ�  �?�Q�  �<¹t  �C�B\  �>���  �C��  �=P�  �@�]  �B�Ǯ  �@XQ�  �A1�  �B�ff  �B�(�  ���&    @*r1�G�  ��<  ��B�H  �?���  �C�e�  �C��  �C�    ���v�  �AQ�  �>z�`AG�  �>���A�R  �<�� <L��  ��< �<  !��<       ���9X  �@�7     B}H\q�� C��U��D� $�=��t  ��=�}  �?ճ� ��^� �?� � �A���  �C��  �A�%D  �B��j  �=l�  �>8��  �9��  �=�F  �?��:  �<�"�  �C���  �>p�Q  �B��  �<���  �����  �B��+  �?���  �A�  �B�33  �B�33  ��t�    @����  ��<  ��	�7  �?��u  �C�u�  �C�pj  �C�Y�  ���R  �A ff  �>p��@�\)  �>�\0AG�  �<#� <�c  ��< �<  !��<       ����-  �@�     B}H\sp  C�� �<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  ��<  �B��  �@
=  �A33  �B���  �B��
  ��	7L    ?n�;�{  ��<  ���\)  �@U�-  �C���  �C�Yd  �C��f  �?9�  �@���  �>���@�^5  �>�Z @�z�  �=8R ;�u�  ��< �<  !��<       �=ȴ9  �@��     B}H\t+� C�«>�}�  ����d  �:�(7  ���-�  �=��  �>�܈ �A���  �C�E  �A���  �B��P  �<��  �=
�  �9���  �=���  �@n�  �<�S7  �Cu�4  �?}xm  �B��h  �;*��  ��<  �B��  �@ b  �A!G�  �B�33  �B���  ���    @D(���\)  ��<  �AiG�  �@֟�  �C�U  �C��i  �C���  �A���  �@��
  �>\(�@�\  �>\ @���  �=8R ;��8  ��< �<  !��<       �@N��  �@��     B}H\u�  C��U>W� $��B  �=+�0  ��U�# �@z�[ �=�� �A��r  �C��  �A��  �B���  �=���  �=��  �;�#  �=�)�  �AU9\  �:�	S  �C��e  �@��{  �B祔  �=���  �A��  �B��+  �@w|�  �AJ�\  �B�ff  �Bsp�  ��/G�    A�Q.�
=  ��<  �B�33  �A���  �C�G�  �C���  �C�ff  �B�33  �AU�  �=���APff  �?w
8AY  �=�\@�<  !��< �<  !��<       �A�{  �@�     B}H\v�� C�� ?�&�  ��i$m  �>F�  ����  �AcdD  �>�q� �A�	  �C��Q  �A��K  �B�±  �=C�  �=���  �<+0�  �>�_  �B���  �<2h(  �C[_�  �@F  �B~�  �=-9�  ��A&�  �B��P  �@�+  �A`��  �B�ff  �Bj�R  �B�33    A�34B=q  ��<  �C(�f  �B33  �C�O(  �C�&2  �C���  �C�@   �A��\  �>�G�A���  �>���A���  �?��=m  ��< �<  !��<       �B�\  �@�     B}H\xR  C�ʫ?@�  �>�C  �>�  �A'�  �AtD|  �?�� �A��  �C�_:  �A�۟  �B�Ï  �?A��  �=�Ϝ  �=Ľ  �?���  �B���  �<1v�  �Co�*  �@+	�  �BoT�  �>'��  ��/�  �B�׎  �@ڟ�  �A���  �B�    �BV�\  �B���    A���Bř�  ��<  �C��   �BS
=  �C��(  �C��G  �C���  �D�3  �A��  �>�Q�A�ff  �>z�A��\  �>���=u׈  ��< �<  !��<       �Bmz�  �@�}     B}H\y� C��U@
c�  �@)v  �>q��  �Bl�  �A��7  �@8K  �A��  �D �"  �A~i�  �B��j  �>��  �>٢�  �=�X�  �>�`  �B�>�  �;��7  �Ch6  �@��  �BPh�  �>�+�  ��ژT  �B��  �Aff  �A��  �B���  �BA�  �CffB   �< C.ff  ��<  �C�Y�  �B�
=  �C���  �C��V  �C�@   �D4@   �A��  �?fpA���  �=�
�A��H  �>W
@=d�t  ��< �<  !��<       �B���  �@�^     B}H\z�  C�� ?���  �@��<  �>Cj  �B�[F  �A��u  �@Nb  �A���  �D5c  �Az9  �Bù�  �>�h�  �>�Z�  �=���  �>�s�  �B�  �;���  �Cd�  �?�'�  �BCA  �>bԉ  ��a�  �B�ɺ  �A%  �A��
  �BvQ�  �B,(�  �BNp�    B�aHC{�  ��<  �C�@   �B�    �C�PA  �C�Z  �C�s3  �Da�   �A�
>  �?L��A��  �>��@A��
  �=�\ =xa�  ��< �<  !��<       �B̙�  �@�?     B}H\{~� C�ҫ?���  �@�h�  �>;6N  �C3�  �A�Q�  �@H�  �A�+�  �D �  �Av�  �BÞ   �?�\  �?Wo�  �=�¢  �?UW�  �Bs�d  �<��  �Co�_  �?�3�  �Bn�|  �>=6Z  ���|�  �Bù�  �AP    �A�
=  �BD��  �B��  �C�fB   �< C�L�  ��<  �D�f  �B���  �C���  �C���  �C�ٚ  �D��   �A�Q�  �?���A�{  �>��Aͅ  �>��=��/  ��< �<  !��<       �B�33  �@�      B}H\}4  C��U?���  �@�p  �>�G  �C��  �Bd{�  �@�o�  �A��  �C���  �Ap�  �BÒ�  �?�[�  �?γN  �>G;�  �?%4�  �BҒ!  �;��  �CuІ  �@
�P  �B��  �>��k  ��B%  �Bí�  �Au�  �A�
=  �A�z�  �A�\  �C�B   �< C��   ��<  �DFf  �B���  �C��  �C��J  �C�&f  �D�@   �A�\)  �?�30A�\  �>�\ A�ff  �?(��=�2�  ��< �<  !��<       �B���  �@�     B}H\}� C�� ?m��  �A�1  �>��6  �CA�|  �BI�T  �@ו�  �A��_  �C��~  �Ao4E  �BÑ  �?���  �?�r*  �>b��  �?"��  �A��  �;6�t  �CpSo  �@j�  �BqM�  �>�!  �� �  �Bå`  �A��  �A���  �A�  �A�
=  �B�ffB   �< Cˀ   ��<  �D*    �B�ff  �C�{�  �C��'  �C��  �D�`   �A��  �?5A�R  �>�p�A��
  �?��=Ƨ�  ��< �<  !��<       �B�33  �@��     B}H\�  C�ګ@b  �A
�  �>�e�  �CN�  �Bc@�  �@��  �A���  �C��  �AmCf  �B�|�  �?�-F  �?�"�  �>x�L  �?%�  �Bx�  �;
[  �CmX  �?�zO  �BdU_  �>叱  ���  �BÇ�  �A���  �A�ff  �A�
=  �A��  �C  B   �< C�&f  ��<  �D5@   �B�33  �C��U  �C��_  �C���  �D�    �A��H  �>z�A�Q�  �>�z�A�G�  �>   =�g8  ��< �<  !��<       �Cff  �@��     B}H\�`� C��U@��  �A*n_  �>�%  �C�$  �BR*  �@�6�  �A��  �C�n�  �Ak�W  �B�V�  �@��  �@4|A  �>���  �?H��  �A��Q  �:�̝  �Ciʆ  �@��  �BW*  �?Tt  ��AA  �B�ff  �A�=q  �A�  �A��  �A���  �B�      A���C�ٚ  ��<  �D;�   �B�ff  �C�f  �C�	�  �C��  �D��   �B (�  �?z�B�  �>ff�B �R  �<�� =���  ��< �<  !��<       �C    �@�     B}H\�  C�� ?��  �A$��  �>�z�  �Cx��  �BN�Y  �@�I�  �A��A  �C��  �Aj��  �B�5p  �?�$P  �?�ǩ  �>�.�  �?Q�  �A�a�  �:��m  �Cg^�  �@@�  �BM{*  �>���  ���  �B�E�  �A�
=  �A�  �A���  �A�Q�  �B��f    A�ffC�Y�  ��<  �D>�   �Bܙ�  �C�  �D   �C�s3  �D�`   �B#�  �?~�@B�\  �?�B �H  �>���=͞�  ��< :��4  ��<       �C�f  �@�B�    B}H\�р C��?�Z  �A5�  �>�L	  �CiC�  �B)��  �@اq  �A�c�  �Cݿ�  �Aj�  �B�+  �?�O`  �?�=  �>c�N  �?Jka  �A���  �:нR  �Ce��  �?��  �BF�j  �>�  ��,�  �B�!H  �A��  �A�33  �A��  �A��\  �Bo    A2Q�C�s3  ��<  �D=�   �B�33  �C��g  �D   �C���  �D�    �Bz�  �?�f`B�  �?�  B�  �>�G�=�;�  ��< 9��E  ��<       �B���  �@�     B}H\��  C��U?ȭ�  �A"  �>�d�  �Cm��  �B+r�  �@�f�  �A�  �C�|  �Ai��  �B��  �?�1  �@;�  �>\�z  �?aF-  �A�h  �:���  �Ce�=  �?��  �BG��  �>̖  ����f  �B��  �A��  �A�  �A�  �A���  �A�\)    ?H��C��  ��<  �D8@   �B���  �C�p�  �D�.  �C��3  �D��   �B��  �?�p�B    �?Y��B�3  �>�� =��  ��< :��  ��<       �B�    �@�#�    B}H\�B� C�� ?�j�  �A,�  �>w�t  �CJh  �BD�  �@��  �B �  �Cޢ]  �Aie4  �B��-  �?̡�  �@e;�  �>Cc�  �?H�P  �B.a�  �:Ґ�  �Ch�I  �@��  �BS�"  �>�Э  ���4  �B��f  �A���  �A��H  �A���  �A���  ��У�    @ߙ�C�ff  ��<  �D/@   �B�33  �C�*�  �D!z  �C�33  �D��   �B�  �?�  B8R  �?��B�  �>(��=Ϫ�  ��< ;a�  ��<       �B�33  �@�     B}H\��  C��?�j  �A�  �>Z��  �COE  �B�  �@�y�  �B �I  �C���  �Ai<?  �B��  �@4z�  �@�  �>gI�  �?,+�  �A��  �:��  �CZ�  �@ �c  �B~<  �>�V�  ����  �B��J  �A�p�  �A�\)  �A��  �A�33  ���Q�    A��C�L�  ��<  �D"S3  �B�    �C��  �D q�  �C��  �D��   �B  �?!G�BL�  �>�f�B�  �>�� =�y�  ��< :�,(  ��<       �B�    �@��    B}H\��� C��U?���  �@޸	  �>�m  �C)�  �A��  �@�u�  �B �)  �C�&Y  �AiKU  �B¯�  �?�?U  �?�`B  �>v��  �?�B  �A�Bx  �:�G  �CZbs  �?��  �B��  �>�-/  ��0��  �B»�  �A��R  �A�Q�  �A��  �A�z�  ��x'    @��oC���  ��<  �D��  �B���  �C�Pw  �C��  �C�    �D��   �BǮ  �?� B�
  �>#� B\  �>�f�=��Z  ��< :*d�  ��<       �B���  �@�u     B}H\�i  C�� ?�B�  �@�oC  �>��  �CT  �A�4�  �@ر]  �Bg  �C�?�  �AiN�  �Bc  �?�f�  �@ d�  �>b�  �>�d  �B��  �:A�  �CV7X  �@v  �B�_  �>ƞ�  ��
�  �B±�  �A��H  �A���  �A��H  �A��  ���A�    A�^C�33  ��<  �C��f  �B���  �C��   �C���  �C���  �Dj@   �B=p  �>��B�  �=�� B�H  �>�R =���  ��< 8ѷ  ��<       �B�    �@��    B}H\�$� C��?��R  �@�.  �=��  �BښB  �Ag�G  �@���  �B h�  �C�{�  �Aic�  �B�  �@�  �?�r�  �>Z�  �>{v  �B�θ  �9�.�  �CO�  �@G  �Aؔ	  �>�L  ��E�  �B¢N  �A��R  �A�\  �A��H  �A�  ��d��    AN��CBL�  ��<  �C�    �B�\  �C��Y  �C��  �C���  �D<�   �B�R  �<�� B�f  �>ff�B�  �>��=��  ��< �<  !��<       �B���  �@�V     B}H\��  C��U?X�  �@+$J  �>c�  �B���  �A���  �@���  �B ��  �C�}�  �Ai%A  �B=  �?���  �?���  �>E2U  �=��  �A�kZ  �9N��  �CO@  �@�:  �A� d  �>�Z�  �®�  �B�  �A��  �A��  �A�=q  �A���  ��p�    @��LB�    ��<  �C�ٚ  �BO��  �C��c  �C�dF  �C�33  �D
@   �A�G�  �>�=�A��
  �>`B �{  �>�) =�  ��< �<  !��<       �BV�R  �@�ƀ    B}H\��� C�� ?e |  �?t��  �=s  �A���  �@��  �@��  �A�
�  �C��  �Ai��  �B�  �?���  �?b'�  �>G�<  �= ^@  �A>'  �9h�  �CF|.  �@  �A��w  �>��=  ��i`  �B7  �A��\  �A�    �A�\)  �A��  �� �    A5ҐA�ff  ��<  �CAL�  �B	�
  �C�l�  �C��  �C�ٚ  �C�ff  �A�\  �>�
@A�z�  �>�� A���  �>�=�=��  ��< �<  !��<       �B�
  �@�7     B}H\�K  C�����  �� fz  ��3�.  ��B�  ����s  �@��  �A���  �C�*�  �AjW�  �B�y�  �?�Ӟ  �?o�  �>%r�  �=���  �C��  �9ĦQ  �C>��  �@t  �A(��  �>��A  �C��F  �B�  �A���  �A���  �B  �Ař�  ��]+     AA\���{  ��<  �B�    �A�    �C�N�  �C��  �C�&f  �C�  �A��  �>���A�  �>�G�A�p�  �>(��=f{_  ��< �<  !��<       �A�(�  �@���    B}H\�� C��U��B   ���=V  ���1B  ��c�  ���'�  �@�Z�  �A��  �D �  �Ak�}  �B�t�  �?Z�|  �>��(  �>C,  �>�  �BL��  �:w��  �C6�  �@	�]  �@�T  �>�b�  �B�7z  �B�|�  �A���  �A��  �B�
  �A�p�  �� �    ?�L�\)  ��<  �A.�R  �@3�  �C�ʖ  �C�ث  �C��  �A�p�  �A�  �=�
�A��H  �=�G�A��  �<#� =��  ��< �<  !��<       �@5p�  �